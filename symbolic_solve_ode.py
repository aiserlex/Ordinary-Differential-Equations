import sympy as s

x = s.Symbol("x")
y = s.Function("y")(x)
y1 = s.diff(y, x)
y2 = s.diff(y, x, x)

ode = s.Eq(y2 + y, 4*s.exp(x))
sol = s.dsolve(ode)
s.pprint(sol)
