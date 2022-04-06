# Solow Dynamics
We solve the Solow Growth Model in Julia using Modeling Toolkit. First we load the packages that we will be using

```Julia
using ModelingToolkit, DifferentialEquations, NonlinearSolve, Plots
```
Then we define our parameteres.

```Julia
@parameters t A α s δ g n
```
and our variables

```Julia
@variables k(t) y(t) c(t)
```
Finally we define an operator for the differentiation w.r.t. time

```Julia
D = Differential(t)
```

Now we write the basic equations for our model. All variables are expressed per effective units of labor.

```Julia
eqs = [ y ~ A*k^α,
        c ~ (1 - s)*y,
        D(k) ~ y - c - (δ + g + n)*k ]
```
We can then tell Julia that this is an ODE

```Julia
@named Solow = ODESystem(eqs)
```
If we try to solve the model now. Julia will try to solve and ADE (Alegabraic Differential Equation) but we know that the Solow model can be expressed as a single differential equation which only depends on *k(t)*. So we can ask Julia to simplify things for us.

```Julia
Solow = structural_simplify(Solow)
```
Notice that before we simplified our model, Julia treated Solow as a model with 3 state variables and 3 equations. After we simplication it became a model with one state variable and one equation.

If we want to simulate how our economy will evolve over time we need to define our parameters

```Julia
p = [A => 1.0, α => 0.5, s => 0.2, δ => 0.065, g => 0.02, n => 0.015]
```
our initial condition

```Julia
k0 = [k => 2.0]
```
and the time interval

```Julia
tspan = (0.0,50.0)
```

that is, we start at *t = 0* and the solution is computed until *t = 50*. The next step is to set the problem in Julia

```Julia
prob = ODEProblem(Solow, k0, tspan, p)
```

and to solve it

```Julia
sol = solve(prob)
```

We can then plot the solution

```Julia
plot(sol, vars=[k,y,c])
```
![](https://github.com/alerodri1976/Solow/blob/main/Solow_1.png)

# Solow Steady State


ss_eqs = [      y ~ A*k^α,
                c ~ (1 - s)*y,
                0 ~ y - c - (δ + g + n)*k ]

@named ss_sys = NonlinearSystem(ss_eqs, [c,y,k], [A,s,α,δ,g,n])
ss_sys = structural_simplify(ss_sys)
ss_prob = NonlinearProblem(ss_sys,[k=>2.0],p)
ss_sol = solve(ss_prob)

k0 = [k => ss_sol[k]]
prob = ODEProblem(Solow, k0, tspan, p)
sol = solve(prob)
plot(sol, vars=[k,y,c])
p2 = [A => 1, α => 0.5, s => 0.3, δ => 0.065, g => 0.02, n => 0.015]
k02 = [k => sol(10.0,idxs=k)]
tspan2 = (10.0,50.0)
prob2 = ODEProblem(Solow, k02, tspan2, p2)
sol2 = solve(prob2)
plot!(sol2, vars=[k,y,c],xlims=tspan)
```
