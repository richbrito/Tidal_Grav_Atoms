# Tidal Love numbers of Gravitational Atoms
This repository contains a Mathematica package that can be used to compute the tidal Love numbers of gravitational atoms. It refers to the work with arXiv number [24XX.XXXXX](https://arxiv.org/abs/24XX.XXXXX), which we ask to cite if you use this code.

The user need only call the function Lovenumber[n,li,mi,l], where {n,li,mi} are the numbers that describe the gravitational atom, while l is the index of the multipolar field that is being tidally induced on the atom. The expression for the Love number is eq. 67 in the work mentioned above, and it is ultimately computed by using the expressions present in Appendix F.2.

We present in the notebook love_scalar.nb a few examples of applications of our code. We compute some simple cases (l=2 for n=0,1 li=mi=0), a case where the Love number is undefined, and a case that requires longer computing time (l=2, n=2,li=mi=2). This computing time is expected to increase a great deal as the indices increase.

The precision of the integrals we compute was set so as to find a middle ground between precision and reasonable computing time, and it can be changed by editing the .wl file.

## Credit

The code is developed and maintained by Ricardo Arana, Richard Brito and Gonçalo Castro. Please, report any issues to

ricardo.arana@tecnico.ulisboa.pt, richard.brito@tecnico.ulisboa.pt and goncalo.e.castro@tecnico.ulisboa.pt

If you make use of these notebooks for your own publications, please cite:
```
@article{XXX}
