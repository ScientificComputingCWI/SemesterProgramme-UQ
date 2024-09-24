# CWI UQ semester programme 

This is the repository for the CWI Research Semester Programme [Uncertainty Quantification for High-Dimensional Problems.](https://www.cwi.nl/en/events/cwi-research-semester-programs/research-semester-programs-in-2024/uncertainty-quantification-for-high-dimensional-problems/).

# Autumn school (7-11 Oct)

Event page: https://www.cwi.nl/en/events/cwi-research-semester-programs/uncertainty-quantification-for-high-dimensional-problems-autumn-school/

## Getting access to Archer2

The lecture on Wednesday (October 9th) includes exercises on submitting ensemble runs on the [ARCHER2 supercomputer](https://www.archer2.ac.uk/), provided by UKRI, EPCC, HPE Cray and the University of Edinburgh. These resources are made available through the [SEAVEA project](https://www.seavea-project.org/) (Software Environment for Actionable & VVUQ-evaluated Exascale Applications).

It is strongly suggested that all course participants aquire a SAFE account and course project accounts in good time well before the 9th of October.

The course project has been created, namely `ta171 - 241009 SEAVEA training day`,  and is ready for use until the 8th of November 2024. All users of this project must have an account on the SAFE, which is EPCC's web-based administration portal, wherein users can request access to projects on any of EPCC's services. To get an account on SAFE, please visit [1].  **You are required to use an institutional email address, e.g. not gmail.com, and also provide your public ssh key**.

To join the course project, hover over the Projects button along the top banner, and click on 'Request access', then enter 'ta171' to locate our course project, then click Apply.

A quick-start guide on how to use Archer2 can be found here [2], with the full user guide here [3].

NB: Archer2 charges per node, so a 4-core job running for an hour will cost 128 core hours (=1CU, where a CU is 1 node hour).   As such, an ill-formed batch script can be expensive, e.g., request 10 hours on multiple nodes and forget a carriage return on the srun command.

The account code/budget to use in batch scripts is "ta171", how to enter this will be made clear in the lecture.

[1] https://safe.epcc.ed.ac.uk/main.jsp

[2] https://docs.archer2.ac.uk/quick-start/quickstart-users/

[3] https://docs.archer2.ac.uk/user-guide/

## Files

* `forward_UQ/`

# Workshop (11-15 Nov)
