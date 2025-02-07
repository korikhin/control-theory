# Attitude Control

This project examines the implementation of attitude control for a small spacecraft of the CubeSat type. Maneuvering is achieved through a reaction wheel assembly in a pyramidal configuration with an inclination angle of $\alpha = 54.73^\circ$.

The control signal is formed based on the current orientation quaternion and angular velocity of the spacecraft:

$$ \tau = -\gamma\mathrm{sign}\left(\mathrm{scal}\left(\bar{\varphi}\circ\lambda\right)\right)\mathrm{vect}\left(\bar{\varphi}\circ\lambda\right) - \beta\omega $$

Where $\lambda$ and $\varphi$ are the current and final orientation quaternions, respectively.

The control torque is distributed among the reaction wheels based on the following expression:

$$ m = A^+\tau $$

Where $A^+$ is the Moore-Penrose inverse of the direction cosine matrix

$$
A = \begin{pmatrix}
  \sin\alpha &          0 & -\sin\alpha &           0 \\
           0 & \sin\alpha &           0 & -\sin\alpha \\
  \cos\alpha & \cos\alpha &  \cos\alpha &  \cos\alpha
\end{pmatrix}
$$

The mathematical model of the process consists of three main components:

- Quaternion Poisson kinematic equation.
- Euler's equations of rotational motion.
- Linearized reaction wheel model in the form of an aperiodic element with damping coefficient $c_w$.

The model is written as a system of vector ODEs as follows:

$$
\begin{cases}
  \dot{\lambda} &= \frac{1}{2}\lambda\circ\omega \\
  J\dot{\omega} + \omega\times J\omega &= \tau \\
  J_w\dot{\Omega} + c_w\Omega &= m
\end{cases}
$$

See the simulation [results](results/attitude_control_simulation.pdf).
