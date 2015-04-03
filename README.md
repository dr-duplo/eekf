# EEKF - Embedded Extended Kalman Filter

This project implements an Extended Kalman Filter in C intended for the use in embedded applications.
The main features are:

- small implementation
- simple C interface using callbacks for state transition and measurement prediction functions
- usable for nonlinear (extended) and linear Kalman Filter cases
- no dynamic memory allocation
- dedicated minimal matrix computation module
- efficient filter computation using Cholesky Factorization
- separated prediction and correction steps
- input and measurment dimension are allowed to change between steps

## What is a Kalman Filter?

With a Kalman Filter one can estimate the internal hidden states of a process/system by measuring only the visible outputs.
This is widely used with Inertia Measurement Units (IMU) to do sensor fusion or dead reconing. As long as no measurement is available the filter will predict the current state of the system. Once a measurement is available it will update the estimate. This is known as prediction and correction step. For more information please refer to http://en.wikipedia.org/wiki/Kalman_filter.

## So what is an Extended Kalman Filter?

In the linear filter case the states from one time step to the next are linear related. That means their exists a constant matrix which expresses that. Also the the output of the system is linear related to the internal states. What if the relations are nonlinear? In this case the relation between the consecutive states and the output is described by nonlinear functions, for which equations of the filter are not solvable. The solution is to linearize the nonlinear functions at the current state and to apply the equations of the linear Kalman. This is what the Extended Kalman filter is all about. For more information please refer to http://en.wikipedia.org/wiki/Kalman_filter.

## What the implementation does

The implementation provides all Kalman Filter computations except for the state prediction function f and the measurment prediction function h. The user has to implemnt these by providing the state and measurement prediction computation and the derivation of the functions with respect to the current filter state (Jacobians). This is done in callbacks. You can use the interface for linear Kalman filter case too. Just let the callbacks return constant Jacobians. The example program shows this approach.




