#src/deeprbp/explainability_module/deeplift_handler.py

import torch
from captum.attr import DeepLift, NeuronDeepLift
import numpy as np
import pandas as pd
from tqdm import tqdm

from ..util.logger import Logger

import warnings
warnings.filterwarnings('ignore')
#warnings.filterwarnings("ignore", message="Input Tensor .* did not already require gradients")
warnings.filterwarnings("ignore", message="Setting forward, backward hooks and attributes on non-linear activations")
    
class DeepLiftHandler:
    """
    DeepLiftHandler class for managing the DeepLIFT explanation process within the DeepRBP explainability module.

    This class prepares reference data, computes attribution scores using the DeepLIFT method,
    and aggregates results to provide insights into the contributions of RNA-binding proteins (RBPs)
    to gene expression.

    Attributes:
        logger (Logger): Logger instance for tracking the progress and results of the explanation process.
        config_explain (ConfigParser): Configuration parser instance for the explanation process, containing settings and paths.
        dataset (DeepRBPExpressionDataset): An instance of DeepRBPExpressionDataset containing RBP, gene, and transcript expression data.
        verbose (int, optional): Verbosity level for logging. Controls the amount of information printed during training.
                                 - 0: No logging (suppress device information and other logs).
                                 - 1: Basic logging (show training progress and essential logs).
                                 - 2: Detailed logging (show device information).
                                 Default is 1.  
    """
    def __init__(self, config_explain, dataset, model, verbose=1):
        """
        Initialize the DeepLiftHandler with the model, dataset, explanation configuration, and gene mappings.

        Parameters:
            model: The trained machine learning model to be explained.
            config_explain (ConfigParser): Configuration parser instance containing settings and paths for the explanation process.
        """
        self.logger = Logger(verbose=verbose)
        self.config_explain = config_explain
        self.dataset = dataset
        self.model = model
        self.model.eval()  # BN en eval para atribuciones estables NEW
        self.logger.log("✅ DeepLiftHandler initialized successfully.", level=1)
 
    def prepare_rbp_reference_tensor(self):
        """
        Prepare reference RBP tensor based on the configuration.
        Returns:
            reference_rbp_tensor (Tensor): Tensor used as a reference for RBP expression.
        """
        # Determine reference tensor based on configuration
        reference_type = self.config_explain.get('reference_data')
        if reference_type == 'median_reference':
            self.logger.log("[prepare_rbp_tensors] Using median of RBP expression as reference.")
            reference_rbp_tensor = torch.tensor(np.median(self.dataset.features['scaled_rbp_df'], axis=0), dtype=torch.float32)
            reference_rbp_tensor = torch.reshape(reference_rbp_tensor, (1, reference_rbp_tensor.size(0)))
            self.logger.log("[prepare_rbp_tensors] Median reference tensor created.")
   
        elif reference_type == 'knockout_reference':
            self.logger.log("[prepare_rbp_tensors] Using knockout (zeroed) RBP expression as reference.")
            reference_rbp_tensor = torch.zeros(1, self.dataset.features['scaled_rbp_df'].shape[1], dtype=torch.float32)
            self.logger.log("[prepare_rbp_tensors] knockout reference tensor created.")
  
        elif reference_type == 'half_reference': # esta no me acaba de convencer, no habria que hacer la mitad de la expresión actual??
            self.logger.log("[prepare_rbp_tensors] Using half-zeroed RBP expression as reference.")
            reference_rbp_tensor = torch.ones(1, self.dataset.features['scaled_rbp_df'].shape[1], dtype=torch.float32) * 0.5
            self.logger.log("[prepare_rbp_tensors] Half-zeroed reference tensor created.")
        else:
            self.logger.error(f"[prepare_rbp_tensors] Unknown reference type: '{reference_type}'", ValueError)
            raise ValueError(f"Unknown reference type: {reference_type}")
        return reference_rbp_tensor

    def compute_attribution_scores(self, reference_rbp_tensor):
        """
        Compute attributions (B,R) for each transcript k, 
        using the YAML 'target_mode': "final" | "logit".

        Parameters:
            reference_rbp_tensor (torch.Tensor): The reference RBP tensor.
            
        Returns:
            list: A list of tensors [(B,R), ...] of length T.
        """
        mode = str(self.config_explain.get('target_mode', 'final')).lower().strip()
        if mode not in ("final", "logit"):
            raise ValueError("config target_mode must be 'final' or 'logit'.");  

        scaled_rbp_tensor = self.dataset.features['scaled_rbp_df'] # (B, R)
        gene_tensor = self.dataset.features['gene_df'] # (B, T)
        T = self.dataset.features['isoform_df'].shape[1]

        # Explainer instance depending on the mode:
        if mode == "final":
            explainer = DeepLift(self.model)
            use_neuron = False
            self.logger.log("[compute_attribution_scores] mode=final ⇒ DeepLift(self.model), target = transcript index in log2(tpm+1).", level=1)
        else:
            # LOGITS (pre-sigmoid): Linear Layer (128->T)
            logits_layer = self.model.abundance_estimator[9]
            explainer = NeuronDeepLift(self.model, logits_layer)
            use_neuron = True
            self.logger.log("[compute_attribution_scores] mode=logit ⇒ NeuronDeepLift(..., Linear(128→T)), neuron_selector = transcript index.", level=1)

        list_batch_scores = []

        for k in tqdm(range(T), desc="Calculating scores", unit="node"):
            if use_neuron:
                # Attribution to the logit of transcript k with respect to INPUTS (RBPs)
                scores = explainer.attribute(
                    inputs=scaled_rbp_tensor,          # (B, R)
                    baselines=reference_rbp_tensor,    # (1, R)
                    neuron_selector=k,                 # 0..T-1
                    additional_forward_args=gene_tensor,
                    attribute_to_neuron_input=False    # neuron output (logit)
                )  # -> (B, R)
            else:
                # Attribution to the final output of the model (post-sigmoid + your mix with gene in forward) (default behaviour)
                scores = explainer.attribute(
                    inputs=scaled_rbp_tensor,          # (B, R)
                    baselines=reference_rbp_tensor,    # (1, R)
                    target=k,                          # 0..T-1
                    additional_forward_args=gene_tensor
                )  # -> (B, R)
            list_batch_scores.append(scores)
        return list_batch_scores

    def calculate_scores_transcript_level(self):
        """
        Calculate attribution scores at the transcript level using the DeepLIFT method.

        Use target_mode from the YAML:
            - final: DeepLift(self.model) with target=k
            - logit: NeuronDeepLift(..., Linear(128->T)) with neuron_selector=k
        Reduce (B) according to batch_reduction_method and return TxRBP.
        Optionally returns df_per_sample if `save_per_sample_scores` is True in the config.
   
        Returns:
            pd.DataFrame: A DataFrame containing the calculated DeepLIFT scores at the transcript level.
        """
        # Prepare RBP reference tensor based on the selected configuration
        self.logger.log("🔧 Preparing RBP reference tensor...", level=1)
        reference_tensor = self.prepare_rbp_reference_tensor()  
        self.logger.log("✅ RBP reference tensor prepared successfully.", level=1)

        # Compute attribution scores for batch
        self.logger.log("🔢 Computing attribution scores for batch...", level=1)
        list_batch_scores = self.compute_attribution_scores(reference_tensor) # K=T, each tensor (B,R)
        self.logger.log(f"✅ Computing attribution scores for batch.", level=1) 

        # Optional: assemble per-sample DataFrame (no I/O)
        df_per_sample = None
        if bool(self.config_explain.get('save_per_sample_scores', False)):
            self.logger.log("🧩 Assembling per-sample attribution DataFrame (TxR×B)...", level=1)
            stack = np.stack([t.detach().cpu().numpy() for t in list_batch_scores], axis=0)  # (T,B,R)
            T, B, R = stack.shape
            arrays = [[f"S{b}" for b in range(B)] * R, list(self.dataset.rbp_names) * B]
            multi_cols = pd.MultiIndex.from_arrays(arrays, names=("Sample", "RBP"))
            df_per_sample = pd.DataFrame(
                stack.reshape(T, B * R),
                index=self.dataset.trans_names,
                columns=multi_cols
            )
            self.logger.log("✅ df_per_sample DataFrame created in memory.", level=1)

        # Reduce the batch dimension to obtain final scores (RBP x T)
        self.logger.log("🔽 Reducing batch dimension to aggregate scores...", level=1)
        # rows = trasncripts; columns = RBPs
        df_scores_TxRBP = self.reduce_batch_dimension(
            list_batch_scores, 
            row_labels=self.dataset.trans_names, 
            col_labels=self.dataset.rbp_names
        ) 
        self.logger.log("✅ Batch dimension reduction complete. Scores aggregated.", level=1)
        return df_scores_TxRBP, df_per_sample
  
    def calculate_scores_hidden_layer(self): # new
        """
        Attributions from RBPs -> neurons of the last hidden layer (after the 128-unit ReLU).
        Uses LayerDeepLift on self.model.abundance_estimator[8] (the final ReLU).
        """
        # objective layer (final ReLU of the 128 block)
        # abundance_estimator = [0]Lin(?,1024), [1]BN, [2]ReLU,
    #                        [3]Lin,       [4]BN, [5]ReLU,
    #                        [6]Lin(1024,128), [7]BN(128), [8]ReLU,
    #                        [9]Lin(128, out), [10]Sigmoid
        layer_module = self.model.abundance_estimator[8]
        # Hidden dimension (output of the preceding Linear layer)
        hidden_dim = self.model.abundance_estimator[6].out_features  # 128
        #
        reference_tensor = self.prepare_rbp_reference_tensor()
        rbp_inputs = self.dataset.features['scaled_rbp_df'] # (B, 1348)
        gene_inputs = self.dataset.features['gene_df'] # not used here, but kept for a consistent signature
        #
        neuron_explainer = NeuronDeepLift(self.model, layer_module)
        list_batch_scores = []
        for neuron_idx in tqdm(range(hidden_dim), desc="Calculating HL scores", unit="neuron"):
            scores = neuron_explainer.attribute(
                inputs=rbp_inputs,              # shape (B, R)
                baselines=reference_tensor,     # shape (1, R)
                neuron_selector=neuron_idx,     # 0..H-1
                additional_forward_args=gene_inputs,
                attribute_to_neuron_input=True      # <--- pre-ReLU (entrada de la ReLU)
            ) # -> Tensor (B, R)
            list_batch_scores.append(scores)
        #
        hl_labels = [f"H{n}" for n in range(hidden_dim)]
        return self.reduce_batch_dimension(
            list_batch_scores, 
            row_labels=hl_labels, # rows = neurons HL
            col_labels=self.dataset.rbp_names # cols = RBPs
        )

    def reduce_batch_dimension(self, list_batch_scores, row_labels, col_labels, max_zero_std_logs: int = 20):
        """
        Reduce batch dimension (B) over a stack (K, B, R) -> (K, R)
        K = #targets (transcripts or HL neurons), R = #RBPs.

        Parameters
        ----------
        list_batch_scores : list[Tensor]
            List of tensors with shape (B, R), one per target (T or HL).
        row_labels : list[str]
            Row labels (len == K) -> trans_names or H0..H{K-1}.
        col_labels : list[str]
            Column labels (len == R) -> rbp_names.
        max_zero_std_logs : int
            Max number of (row, col) pairs with std == 0 to log (avoid flooding).

        Returns
        ----------
        pd.DataFrame (K x R)
        """
        method = self.config_explain.get('batch_reduction_method')
        if method not in ['t-statistic', 'sum_scores']:
            raise ValueError(f"Invalid batch reduction method: {method}")
            
        # (K, B, R)
        stack = np.stack([t.detach().cpu().numpy() for t in list_batch_scores], axis=0)
        print("stack.shape:", stack.shape)

        # ===================== DEBUG EJEMPLO CONCRETO =====================
        debug_row = "ENST00000610897"
        debug_col = "ENSG00000130764"

        if (debug_row in row_labels) and (debug_col in col_labels):
            r_idx = row_labels.index(debug_row)
            c_idx = col_labels.index(debug_col)

            v = stack[r_idx, :, c_idx]   # vector (B,) con los scores por batch
            print(f"[DEBUG] Raw batch scores for ({debug_row}, {debug_col}): {v}")

            if method == "t-statistic":
                m = v.mean()
                s = v.std()
                n = v.shape[0]
                if s != 0:
                    t_debug = m / (s / np.sqrt(n))
                else:
                    t_debug = 0.0
                print(f"[DEBUG] mean={m}, std={s}, n={n}, t-stat={t_debug}")
            else:  # sum_scores
                print(f"[DEBUG] sum_scores={v.sum()}")
        else:
            print("[DEBUG] Pair (ENST00000610897, ENSG00000130764) not found in labels.")
        # =================================================================

        if method == 't-statistic':
            mean = np.mean(stack, axis=1)  # (K, R)
            std  = np.std(stack, axis=1)   # (K, R)
            n    = stack.shape[1]
            with np.errstate(divide='ignore', invalid='ignore'):
                result = np.where(std != 0, mean / (std / np.sqrt(n)), 0.0)  # (K, R)
            
            # ── NEW: zero-std diagnostics, generic for TxRBP or HL×RBP
            zero_std_mask = (std == 0)
            num_zero = int(zero_std_mask.sum())
            
            if num_zero > 0:
                self.logger.log(f"[reduce_batch_dimension] Found {num_zero} positions with std == 0.", level=1)
                idxs = np.argwhere(zero_std_mask)
                for i, (r, c) in enumerate(idxs[:max_zero_std_logs]):
                    self.logger.log(f"[reduce_batch_dimension] Zero std at Row: {row_labels[int(r)]}, Col: {col_labels[int(c)]}.", level=1)
                if num_zero > max_zero_std_logs:
                    self.logger.log(f"[reduce_batch_dimension] ... and {num_zero - max_zero_std_logs} more.", level=1)
       
        else:  # sum_scores
            result = np.sum(stack, axis=1)  # (K, R)

        K, R = result.shape
        if len(row_labels) != K:
            raise ValueError(f"[reduce] Row labels length {len(row_labels)} != K {K}")
        if len(col_labels) != R:
            raise ValueError(f"[reduce] Col labels length {len(col_labels)} != R {R}")

        # =============== DEBUG: VALOR FINAL COLAPSADO =====================
        if (debug_row in row_labels) and (debug_col in col_labels):
            r_idx = row_labels.index(debug_row)
            c_idx = col_labels.index(debug_col)
            print(f"[DEBUG] Collapsed score result[{debug_row}, {debug_col}] = {result[r_idx, c_idx]}")
        # =================================================================
        return pd.DataFrame(result, index=row_labels, columns=col_labels)
