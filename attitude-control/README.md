# Attitude Control

This project examines attitude control implementation for a CubeSat-type spacecraft. Maneuvering is achieved through a pyramidal reaction wheel assembly with a skew angle of $\alpha = 54.73^\circ$.

## Control Law

The control signal is formed from the spacecraft's current orientation quaternion and angular velocity:

$$ \tau = -\gamma\mathrm{sign}\left(\mathrm{scal}\left(\bar{\varphi}\circ\lambda\right)\right)\mathrm{vect}\left(\bar{\varphi}\circ\lambda\right) - \beta\omega $$

where $\lambda$ and $\varphi$ are the current and target orientation quaternions, respectively.

## Reaction Wheel Distribution

The control torque is distributed among the reaction wheels using:

$$ m = A^+\tau $$

where $A^+$ is the pseudo-inverse of the wheel allocation matrix:

$$
A = \begin{pmatrix}
  \sin\alpha &          0 & -\sin\alpha &           0 \\
           0 & \sin\alpha &           0 & -\sin\alpha \\
  \cos\alpha & \cos\alpha &  \cos\alpha &  \cos\alpha
\end{pmatrix}
$$

## Mathematical Model

The system dynamics are described by the following system of ODEs:

$$
\begin{cases}
  \dot{\lambda} &= \frac{1}{2}\lambda\circ\omega \\
  J\dot{\omega} + \omega\times J\omega &= \tau \\
  J_w\dot{\Omega} + c_w\Omega &= m
\end{cases}
$$

representing the quaternion Poisson kinematic equation, Euler's rotational dynamics, and the linearized reaction wheel model (an aperiodic element with damping coefficient $c_w$), respectively.

See the simulation [results](results/attitude_control_simulation.pdf).
