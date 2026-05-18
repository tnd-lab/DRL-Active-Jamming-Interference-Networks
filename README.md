# Deep Reinforcement Learning-Enabled Anti-Reactive Jamming Strategies for Mobile Cell-Free Interference Networks  
**Tri Nhu Do**, **Georges Kaddoum**, and **François Leduc-Primeau**  
_IEEE Canadian Conference on Electrical and Computer Engineering (CCECE)_  
Canada, May 2026  

## Abstract
    In this paper, we propose a new concept of a mobile cell-free interference network and address the collision-free frequency allocation (CFFA) problem in this network under reactive jamming. The CFFA is formulated as a mixed-integer optimization problem with finite state and action spaces, which can be solved using exhaustive search, assuming the availability of the jammer’s channel state information (CSI). Considering a more realistic scenario in which the jammer’s CSI is unknown prior to the attack, we reformulate the CFFA as an interactive Markov decision process and propose a double deep Q-learning (DDQL) approach to minimize the loss caused by the jammer and achieve CFFA for the remaining active transmissions. Within this framework, the network operates as a cognitive radio system that attains real-time situational awareness through DDQL interactions, enabling adaptive responses to dynamic jamming behaviors. By considering reference point group mobility for legitimate users and random walk mobility for the jammer, we demonstrate that the proposed DDQL achieves a near-optimal policy compared to the optimal exhaustive search. We conclude that, in the presence of unavoidable active jamming, there exists a policy for multiple concurrent cell-free transmission pairs that can minimize spectral efficiency loss in a mobile environment.

## Paper
- [pdf version]()
