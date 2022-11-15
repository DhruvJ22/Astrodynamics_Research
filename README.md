# Astrodynamics Research

This repository contains highlights of the Astrodynamics research that I am doing at Purdue University as a member of the Multi-Body Dynamics Research Group (https://engineering.purdue.edu/people/kathleen.howell.1/index.html). 

## Purpose
I have been working on simulating and investigating the CR3BP model since the beginning of 2021. After a painstaking amount of work, I have been able to successfully simulate the model and target known structures in it. After a lot of refactoring sessions I have converted my initial code that I wrote using "Design after the fact" to be closer to "Design for contribution". 

_**This repo serves as an example of my experience with implementation of numerical methods, numerical optimization, numerical linear algebra, production quality coding, and software design principles.**_

## Setup
The above set of programs provides a robust setup for people to play with and understand periodic orbits and the general dynamical flows in CR3BP.

The code is structured smartly to exploit the principles of Object Oriented Programming in the below manner:
![Dhruv Jain CR3BP dynamics](https://user-images.githubusercontent.com/33181026/159374925-6fe2bf57-9155-48f8-9777-7d9618de9e03.png)

Production quality software  design practices were used to design the programs. A list of some of the princples exploited to make the code reliable (this is not a comprehensive list of all the principles used):
1. Structural Patterns
2. Clear naming convention
3. Elaborate comments
4. Single purpose functions
5. Avoided cryptic nesting
6. Functions with similar argument structure 
7. Single Inhertiance and Encapsulation
8. Modular structure - will enable smooth inclusion of Multipls shooter algorithm

The current work requires the widely used numpy, scipy, plotly and matplotlib libraries. 

Checkout my in-progress contributions to [Poliastro](https://docs.poliastro.space/en/stable/)!

## Examples

Any hamiltonion dynamical model usually has 3 kinds of solutions: equilibrium solutions, periodic solutions and quasi-periodic solutions. In addition, these solutions usually exist as families of solutions.

The current set of programs can compute the equillibrium solutions (Libration Points) and periodic solutions (Periodic Orbits) of CR3BP. Furthermore, it is capable of computing families of periodic solutions (Family of Periodic Orbits).

There are an infinite number of periodic orbits but the three general type of families of known symmetry are Halo family, Axial family, Vertical family. 

1. **Halo family orbits** show XZ plane symmetry. Below is an example of Halo family in Earth-Moon system around L3 libration point. _Some of the members were computed using Natural Parameter Continuation in 'z' using XZ plane symmetry with line search for step size. Other members were computed using Pseudo-arc length continuation using XZ plane symmetry with line search for step size._

![EM_L3_halo](https://user-images.githubusercontent.com/33181026/159376462-be6147ef-a36a-4906-adf0-e8ef35c85bc1.png)

2. **Axial family orbits** show X-axis symmetry. Below is an example of Axial family in Earth-Moon system around L4 libration point. _Some of the members were computed using Natural Parameter Continuation in 'z' by targeting Periodicity with line search for step size. Other members were computed using Pseudo-arc length continuation by targeting Periodicity with pahse constraint and line search for step size._

![EM_L4_axial](https://user-images.githubusercontent.com/33181026/159376511-e3c40cf1-e7b9-47a4-b303-4098de09fbf0.png)

3. **Vertical family orbits** show XZ plane and X-axis symmetry. Below is an example of Vertical family in Earth-Moon system around L2 libration point. _Some of the members were computed using Natural Parameter Continuation in 'x' by targeting XZ plane symmetry with line search for step size. Some of the other members were computed using Natural Parameter Continuation in 'x' by targeting Periodicity with line search for step size. Rest of the members were computed using Natural Parameter Continuation in 'jc' by targeting XZ plane symmetry with line search for step size._

![EM_L2_vertical](https://user-images.githubusercontent.com/33181026/159376601-ca10cf77-0685-46d9-81b8-960a53461c56.png)

The repo consists of multiple other examples to showcase the robustness of the setup and its ability to compute periodic orbits around other libration points and other types of periodic orbits like lyapunov orbits and short period orbits. The examples also serve as "tests" for the various functions and classes.

## Research Focus
The focus of my research is Quasi-Periodic Orbits in the Circular Restricted Three-Body Problem model. The research has many exciting applications, one of the interesting ones being to be able to leverage these dynamical structuresto design host orbits for logistics depots in the cislunar space!

## References 
**Software design principles:**
1. Gamma, E., Helm, R., &amp; Johnson, R. (1998). Design patterns elements of Reusable Object Oriented Software. Addison Wesley. 
2. Wilson, G. (2022). Twelve quick tips for software design. PLOS Computational Biology, 18(2). https://doi.org/10.1371/journal.pcbi.1009809 
3. Silen, P. (2020, March 11). Useful tips for naming your variables. Medium. Retrieved March 21, 2022, from https://betterprogramming.pub/useful-tips-for-naming-your-variables-8139cc8d44b5 
4. 4. Patrick J Mineault & The Good Research Code Handbook Community (2021). The Good Research Code Handbook. Zenodo. doi:10.5281/zenodo.5796873
