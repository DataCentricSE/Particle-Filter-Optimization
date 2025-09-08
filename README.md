# Particle-Filter-Optimization
PFO algorithm for robust optimization. The folder contains the codes for tests on two benchmark functions (2D Branin function and a 1D user-defined Gaussian mixture function) and a practical example of a styrene reactor where profit maximization is the aim.

These codes gives the basis of the work documented in 
Kenyeres, É., Kummer, A., & Abonyi, J. (2025). Improvements of particle filter optimization algorithm for robust optimization under different types of uncertainties. Heliyon, 11(1).
<a href="https://doi.org/10.1016/j.heliyon.2024.e41573">DOI: 10.1016/j.heliyon.2024.e41573</a>

This paper introduces a methodology for handling different types of uncertainties during robust optimization. In real-world industrial optimization problems, many types of uncertainties emerge, e.g., inaccurate setting of control variables, and the parameters of the system model are usually not known precisely. For these reasons, the global optimum considering the nominal values of the parameters may not give the best performance in practice. This paper presents a widely usable sampling-based methodology by improving the Particle Filter Optimization (PFO) algorithm. Case studies on benchmark functions and even on a practical example of a styrene reactor are introduced to verify the applicability of the proposed method on finding robust optimum, and show how the users can tune this algorithm according to their requirement. The results verify that the proposed method is able to find robust optimums efficiently under parameter and decision variable uncertainties, as well.

The folder contains the following codes:
<ul>
<li>pfo.m: general PFO for global optimization</li>
</ul>

Styrene reactor: 
<ul>
<li>main_pfo_styrene_v3.m: experiments on the Robust optimization problem of a styrene reactor</li>
<li>pfo_robust_v1.m: PFO with only parameter uncertainty</li>
<li>untitled10.m: figures</li>
<li>pfo_robust_v2.m: PFO with parameter and decision variable uncertainty (main_styrene_pfo_v3.m/2024.05.09.)</li>
</ul>

Branin function:
<ul>
<li>main_styrene_v4.m: tests on the Branin function and an univariate Gaussian function (some results in sensitivity_dx&dw1.mat)</li>
<li>pfo_robust_v2.m: see above</li>
<li>branin_pso.m: optimization tests on the Branin function by Particle Swarm Optimization (PSO)</li>
</ul>

Gaussian bimodal 1D:
<ul>
<li>main_styrene_v4.m (main_styrene_v4.m/gaussian_bimodal_1d univariate Ns sensitivity, 07.08. --> sensitivity) </li>
<li>pfo_robust_v2.m: see above</li>
<li>benchmark_gaussian_bimodal_1d.m: tuneable multimodal Gaussian function</li>
<li>untitled12.m: expected value curves of benchmark_gauss_bimodal_1d() and violinplot for evaluating results (Violinplot-Matlab-master folder is needed)</li>
</ul>

Results:
<ul>
<li>sensitivity_dx&dw1.mat: some results on the Gaussian function</li>
</ul>
