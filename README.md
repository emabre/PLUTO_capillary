# Copyright
The software contained in this folder and its subfolders (recursively) is licenced under the GNU GENERAL PUBLIC LICENSE (Version 2, June 1991), which may be consulted by reading the file COPYING.
This software is a modification of the code PLUTO v. 4.2 (August 2015), developed by Andrea Mignone, with the contribution of C.Zanni, B.Vaidya, T.Matsakos, G.Musicianisi, P.Tzeferakos, O.Tesilanu (see http://plutocode.ph.unito.it/).
# What is it?
This is a folder containing a slightly modified version of PLUTO. It is meant to work with the other source code files contained inside my project _capillary_, which is also a modifcation of some of PLUTO's source code (and configuration files). In the process of making the present software I was helped by Alberto Marocchino, Stefano Atzeni, Andrea Mignone (I sincerely thank you all).
Everything is provided with no warranty at all. If you use it you do it at your own risk.
# Main changes
I changed the way temperature is computed for thermal conduction (more general than simply p/rho).
I added some code to handle an internal corner in internal boundary.
I made some code to integrate PLUTO with an ADI scheme for thermal conduction and resistivity.
I applied some other minor changes and added some comments.
# Original Readme file
The original Readme file is now named README_original

Emanuele Brentegani
