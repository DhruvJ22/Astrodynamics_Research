# Astrodynamics Research

This repository contains a few highlights of the astrodynamics research that I conducted during my master's at Purdue University between 2021 and 2023 as a member of the [Multi-Body Dynamics Research Group](https://engineering.purdue.edu/people/kathleen.howell.1/index.html). My research primarily focused on innovating frameworks to design transfer trajectories in a multi-body regime between periodic orbits by leveraging quasi-periodic orbits. The programs in this repository model the behavior of a small body in a three-body environment by employing the Circular Restricted Three Body Problem (CR3BP) equations of motion and are capable of computing periodic dynamical structures in the system.

## Purpose
The codebase can be leveraged to simulate the behavior of a small body in a multi-body regime to better understand the general dynamical flows as well as periodic orbits in CR3BP. After multiple refactoring sessions, I have converted the initial code that I wrote using "Design after the fact" to be closer to "Design for contribution".

_**This repo serves as an example of my experience with the implementation of numerical methods, numerical optimization, numerical linear algebra, production quality coding, and software design principles.**_

## Setup
The codebase is structured smartly to exploit the principles of Object Oriented Programming in the below manner:
![Dhruv Jain CR3BP dynamics](https://user-images.githubusercontent.com/33181026/159374925-6fe2bf57-9155-48f8-9777-7d9618de9e03.png)

Production quality software design practices were used to design the programs. A list of some of the principles exploited to make the code reliable:
1. Structural Patterns
2. Clear naming convention
3. Elaborate comments
4. Single purpose functions
5. Avoided cryptic nesting
6. Functions with similar argument structure
7. Single Inheritance and Encapsulation
8. Modular structure - will enable smooth inclusion of Multiple shooter algorithm

The work heavily relies on the following libraries: NumPy, SciPy, Plotly and Matplotlib.

## Examples
The CR3BP model, which is a Hamiltonian dynamical model, commonly possesses 3 types of particular solutions: equilibrium solutions, periodic solutions, and quasi-periodic solutions. In addition, these solutions exist as families of solutions.

The current set of programs can compute:
- Equilibrium solutions (Libration Points) of CR3BP
- Periodic solutions (Periodic Orbits) and their families in CR3BP

Additionally, the CR3BP model possesses two distinct symmetries, and they are leveraged in the codebase to reduce the computational complexity of constructing periodic orbits. The following orbit families exploit the distinct symmetries of the model, and showcase the features of the programs:

1. **Halo family orbits** are governed by XZ plane symmetry. Below is an example of a halo orbit family in the Earth-Moon system around L<sub>3</sub> Lagrange point. Some of the members were computed using _natural parameter continuation_ in 'z' using _XZ plane symmetry_ with _line search_ for step size. Other members were computed using _pseudo-arc length continuation_ using _XZ plane symmetry_ with _line search_ for step size.

![EM_L3_halo](https://user-images.githubusercontent.com/33181026/159376462-be6147ef-a36a-4906-adf0-e8ef35c85bc1.png)

2. **Axial family orbits** are governed by X-axis symmetry. Below is an example of an axial family in the Earth-Moon system around L<sub>4</sub> Lagrange point. Some of the members were computed using _natural parameter continuation_ in 'z' by targeting _periodicity_ with _line search_ for step size. Other members were computed using _pseudo-arc length continuation_ by targeting _periodicity with phase constraint_ and _line search_ for step size.

![EM_L4_axial](https://user-images.githubusercontent.com/33181026/159376511-e3c40cf1-e7b9-47a4-b303-4098de09fbf0.png)

3. **Vertical family orbits** are governed by XZ plane and X-axis symmetry. Below is an example of a vertical family in the Earth-Moon system around L<sub>2</sub> Lagrange point. Some of the members were computed using _natural parameter continuation_ in 'x' by targeting _XZ plane symmetry_ with _line search_ for step size. Some of the other members were computed using _natural parameter continuation_ in 'x' by targeting _periodicity_ with _line search_ for step size. The rest of the members were computed using _natural parameter continuation_ in 'jc' by targeting _XZ plane symmetry_ with _line search_ for step size.

![EM_L2_vertical](https://user-images.githubusercontent.com/33181026/159376601-ca10cf77-0685-46d9-81b8-960a53461c56.png)

The repo consists of multiple other examples to showcase the robustness of the setup and its ability to compute periodic orbits around other Lagrange points and other types of periodic orbits like Lyapunov orbits and short-period orbits.

## References
**Software design principles:**
1. Gamma, E., Helm, R., &amp; Johnson, R. (1998). Design patterns elements of Reusable Object Oriented Software. Addison Wesley.
2. Wilson, G. (2022). Twelve quick tips for software design. PLOS Computational Biology, 18(2). https://doi.org/10.1371/journal.pcbi.1009809
3. Silen, P. (2020, March 11). Useful tips for naming your variables. Medium. Retrieved March 21, 2022, from https://betterprogramming.pub/useful-tips-for-naming-your-variables-8139cc8d44b5
4. Patrick J Mineault & The Good Research Code Handbook Community (2021). The Good Research Code Handbook. Zenodo. doi:10.5281/zenodo.5796873