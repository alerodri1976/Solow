# Solow Growth Model

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

## Simulation

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

This implies that we start at *t = 0* and the solution is computed until *t = 50*. The next step is to set the problem in Julia

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

## Steady State

No we will use a nonlinear solver to find the steady state of the Solow growth model. This model is simple enough so that the steady state can be calculated with pen and paper but let´s solve it numerically. First we define our steady state equations by setting *D(x) = 0*

```Julia
ss_eqs = [      y ~ A*k^α,
                c ~ (1 - s)*y,
                0 ~ y - c - (δ + g + n)*k ]
```
We then tell Julia that this is a nonlinear system

```Julia
@named ss_sys = NonlinearSystem(ss_eqs, [c,y,k], [A,s,α,δ,g,n])
```
and we simplify it

```Julia
ss_sys = structural_simplify(ss_sys)
```

Then we set the problem using our parameters in *p* and we define an initial condition for our solver which is *k = 2*

```Julia
ss_prob = NonlinearProblem(ss_sys,[k=>2.0],p)
```

Finaly we tell Julia to the solution to our problem

```Julia
ss_sol = solve(ss_prob)
```

We can see the steady state values for *y, c* and *k* using symbolic indexing

```Julia
ss_sol[y]
ss_sol[c]
ss_sol[k]
```
## Comparative Statics

Now we show how the solution changes when one of the parameters changes. In this example we will increase the savings rate *s* from 20% to 30% in *t = 10* assuming that we start in a steady state. To do so we will simulate our model under two sets of initial conditions and parameters. First we simulate our baseline model which starts at the steady state.

```Julia
k0 = [k => ss_sol[k]]
prob = ODEProblem(Solow, k0, tspan, p)
sol = solve(prob)
```

Now we setup our alternative world by changing the parameters of the model

```Julia
p2 = [A => 1, α => 0.5, s => 0.3, δ => 0.065, g => 0.02, n => 0.015]
```

the initial condition

```Julia
k02 = [k => sol(10.0,idxs=k)]
```

and the time span

```Julia
tspan2 = (10.0,50.0)
```

We put all this together into a new problem and solve it

```Julia
prob2 = ODEProblem(Solow, k02, tspan2, p2)
sol2 = solve(prob2)
```

To plot our solutions first we plot the solution from the baseline simulation

```Julia
plot(sol, vars=[k,y,c])
```

and the we tell Julia to add the second simulation to the same plot by using the *plot!* command instead of the *plot* function.

```Julia
plot!(sol2, vars=[k,y,c],xlims=tspan)
```
In this case we added the *xlims=tspan* option so that the x-axis range is kept from *t = 0* to *t = 50*

![](https://github.com/alerodri1976/Solow/blob/main/Solow_2.png)

With a little bit of extra commands we can get the labels of the graph in order but that is beyond the scope of this tutorial.
