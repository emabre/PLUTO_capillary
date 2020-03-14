# Copyright
The software contained in this folder and its sub-folders (recursively) is licensed under the GNU GENERAL PUBLIC LICENSE (Version 2, June 1991), which may be consulted by reading the file COPYING.
This software is a modified version of the code PLUTO v. 4.2 (August 2015). PLUTO is developed by Andrea Mignone, with the contribution of C.Zanni, B.Vaidya, T.Matsakos, G.Musicianisi, P.Tzeferakos, O.Tesilanu (see http://plutocode.ph.unito.it/).
# What is it?
This is a slightly modified version of PLUTO, and I used it to produce the numerical simulations of hydrogen-filled capillary discharges that I included in my PhD thesis (http://hdl.handle.net/11573/1361281).
![](ne_evolution_withI.gif)
This software is meant to work together with the source code files contained inside another repository of mine, _capillary_, which is also a modification and an integration of some of PLUTO's source code (and configuration files). In the process of working on the present software I was helped by Alberto Marocchino, Stefano Atzeni, Andrea Mignone (I sincerely thank you all).
Everything is provided with no warranty at all. If you use it you do it at your own risk.
# Main changes
I changed the way temperature is computed for thermal conduction (more general than simply p/rho).
I added some code to handle an internal corner in internal boundary.
I implemented and integrated into PLUTO an _alternating direction implicit_ scheme for thermal conduction and resistivity.
I applied some other minor changes and added some comments.
# Original Readme file
The original Readme file is now named README_original

Emanuele Brentegani
