# **Alec Training Notes** — Physics-Informed Neural Network (PINN)


# 1. Conceptual Overview

A **Physics-Informed Neural Network (PINN)** is a neural network trained using:
1. Observed data
2. A governing physical law (typically a differential equation)

Instead of only fitting data, the model is constrained to satisfy physics everywhere in the domain.

## Governing Equation

We are solving the first-order ODE:

```math
\frac{du}{dx} + u = 0
```

True solution:

```math
u(x) = e^{-x}
```

## Core Idea

We approximate the solution using a neural network:

```math
u_\theta(x)
```

We then enforce:
- Agreement with known data points
- Satisfaction of the differential equation

## Loss Function

```math
\mathcal{L} = \mathcal{L}_{data} + \lambda \mathcal{L}_{physics}
```
Where $\lambda$ is a hyperparameter that balances the two loss components. And where $\mathcal{L}_{data}$ and $\mathcal{L}_{physics}$ represent the data loss and physics loss, respectively. The data loss means the mean squared error between the predicted and true values at the data points, while the physics loss measures how well the predicted solution satisfies the governing differential equation at the collocation points.

### Data Loss

```math
\mathcal{L}_{data} = \frac{1}{N} \sum (u_\theta(x_i) - u(x_i))^2
```
In words, this is the mean squared error between the predicted values from the neural network and the true values at the observed data points.

### Physics Loss

```math
\mathcal{L}_{physics} = \frac{1}{N} \sum \left(\frac{du_\theta}{dx} + u_\theta \right)^2
```
This term penalizes deviations from the governing differential equation at the collocation points, ensuring that the neural network's predictions adhere to the underlying physics.


# 2. Imports

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import torch
    import torch.nn as nn

- `torch`: tensor library + automatic differentiation
- `torch.nn`: neural network modules

Torch is the core library for building and training neural networks, while `torch.nn` provides pre-built layers and activation functions that simplify the construction of neural network architectures. The other imports (`pandas`, `numpy`, and `matplotlib`) are commonly used for data manipulation and visualization, but the core functionality of the PINN relies on PyTorch.


# 3. Input Domains

    x_data = torch.linspace(0,1,5).view(-1,1)
    x_col  = torch.rand(100,1)
    x_test = torch.linspace(0,1,20).view(-1,1)

In the grand scheme of neural network training, these represent the different sets of input points that we will use for training and evaluation:


## Meaning

- `x_data`: the sparse set of points where we have observed data and will compute the data loss.
- `x_col`: the collocation points where we will enforce the physics by computing the physics loss.
- `x_test`: a dense set of points where we will evaluate the trained model to visualize its performance against the true solution.

## `.view(-1,1)`

Reshapes tensor into:

```math
(N,1)
```
Meaning it has `N` rows and 1 column. This is required because `nn.Linear` expects 2D input.


# 4. Target Function

    u_data = torch.exp(-x_data)

Represents:

```math
u(x) = e^{-x}
```


# 5. Neural Network Model

    class PINN(nn.Module):
        def __init__(self):
            super().__init__()

## What `super().__init__()` does

- Initializes PyTorch internal model structure  
- Registers parameters for training  
- Enables gradient tracking  

## Architecture

    self.net = nn.Sequential(
        nn.Linear(1, 32),
        nn.Tanh(),
        nn.Linear(32, 32),
        nn.Tanh(),
        nn.Linear(32, 1)
    )

### Linear Layer

```math
y = Wx + b
```
`nn.Linear` implements a fully connected layer that applies a linear transformation to the input data, where `W` is the weight matrix and `b` is the bias vector. The first layer transforms the 1-dimensional input into a 32-dimensional hidden representation, and the final layer maps the 32-dimensional hidden representation back to a single output.

### Tanh Activation

```math
\tanh(x) = \frac{e^x - e^{-x}}{e^x + e^{-x}}
```

---

## Forward Pass

    def forward(self, x):
        return self.net(x)

Defines:

```math
u_\theta(x)
```


# 6. Optimizer (Adam)

    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

## What Adam Does

```math
\theta_{t+1} = \theta_t - \alpha \frac{\hat{m}_t}{\sqrt{\hat{v}_t} + \epsilon}
```

### Internal Variables

```math
m_t = \beta_1 m_{t-1} + (1-\beta_1) g_t
```

```math
v_t = \beta_2 v_{t-1} + (1-\beta_2) g_t^2
```

```math
\hat{m}_t = \frac{m_t}{1-\beta_1^t}, \quad
\hat{v}_t = \frac{v_t}{1-\beta_2^t}
```


# 7. Training Loop

    for epoch in range(epochs):



## Zero Gradients

    optimizer.zero_grad()

PyTorch accumulates gradients:

```math
\nabla = \nabla + \nabla_{new}
```


## Data Loss

    u_pred_data = model(x_data)
    data_loss = torch.mean((u_pred_data - u_data) ** 2)

```math
\mathcal{L}_{data} = \frac{1}{N} \sum (u_\theta - u)^2
```


## Enable Gradients on Inputs

    x_col_epoch = x_col.clone().detach().requires_grad_(True)

Needed to compute:

```math
\frac{du}{dx}
```


## Forward Pass (Physics)

    u_col = model(x_col_epoch)


## Automatic Differentiation

    du_dx = torch.autograd.grad(
        outputs=u_col,
        inputs=x_col_epoch,
        grad_outputs=torch.ones_like(u_col),
        create_graph=True
    )[0]

Computes:

```math
\frac{\partial u}{\partial x}
```

## Physics Loss

    physics_loss = torch.mean((du_dx + u_col) ** 2)

```math
r(x) = \frac{du}{dx} + u
```

```math
\mathcal{L}_{physics} = \|r(x)\|^2
```


## Total Loss

    loss = data_loss + lambda_phys * physics_loss


## Backpropagation

    loss.backward()

Computes:

```math
\frac{\partial \mathcal{L}}{\partial \theta}
```


## Parameter Update

    optimizer.step()

Applies Adam update rule.


# 8. Evaluation

    u_pred = model(x_test).detach().numpy()
    u_true = torch.exp(-x_test).detach().numpy()

## `.detach()`

- Removes tensor from computation graph  
- Required before converting to NumPy  


# 9. Final Insight

```math
\text{Differential Equation} \rightarrow \text{Optimization Problem}
```

- Neural network approximates solution  
- Autograd enforces physics  
- Loss combines data + PDE  


# 10. Critical Notes

- `requires_grad=True` is essential  
- `create_graph=True` enables higher derivatives  
- PINNs are sensitive to:
  - learning rate  
  - loss weighting  
  - network architecture  