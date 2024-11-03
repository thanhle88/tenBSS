# Tensor-based blind source separation

We proposed novel tensor-based BSS methods, namely TenSOFO and TCBSS. The former is designed for a joint
individual differences in scaling (INDSCAL) decomposition, addressing instantaneous (linear) BSS tasks; while
the latter efficiently performs a constrained block term decomposition (BTD), aligning with the design of
convolutive BSS. 


## Demo
Please run 

+ `demo_joint_decomposition.m`: To see the performance of the proposed joint INDSCAL decomposition of two tensors
+ `demo_type2_BTD.m`: To see the performance of the proposed type2-BTD decomposition
+ `demo_source_separation.m`: To illustrate the performance of tensor-based source seperation algorithms
+ `demo_EMG.m`: To see the performance of our algorithm for EMG source separation
## Reference

This code is free and open source for research purposes. If you use this code, please acknowledge the following paper.

[1] L.T. Thanh, K. Abed-Meraim, P. Ravier, O. Buttelli, A. Holobar. "[*Tensorial Convolutive Blind Source Separation*](https://ieeexplore.ieee.org/document/10447269)". **Proc. 49th IEEE ICASSP**, 2024. 

[2] L.T. Thanh, K. Abed-Meraim, P. Ravier, O. Buttelli, A. Holobar. "[*Joint INDSCAL Decomposition Meets Blind Source Separation*](https://ieeexplore.ieee.org/document/10447387)". **Proc. 49th IEEE ICASSP**, 2024. 
