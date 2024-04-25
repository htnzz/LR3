import sympy as sy


def func(x):
    return x ** 2 - sy.cos(sy.pi * x)


def values(a, b, step):
    table = []
    x = a
    while x <= b:
        table.append((x, func(x)))
        x += step
    return table


def lagrange_polynomial(points, x):
    l = 0
    for i, (x_i, y_i) in enumerate(points):
        l_i = 1
        for j, (x_j, _) in enumerate(points):
            if i != j:
                l_i *= (x - x_j) / (x_i - x_j)
        l += y_i * l_i
    return l


def take_diff(func, x, n):
    new_func = func
    for _ in range(n):
        new_func = sy.diff(new_func, x)
    return new_func


def omega(a, b, step, x):
    res = 1
    while round(a, 2) <= b:
        res *= (x - a)
        a += step
    return res


def main():
    x = sy.symbols('x')
    k = 1
    m = 0
    n = 5
    a = 0.1
    b = 0.6
    step = (b - a) / 3

    points = values(a, b, step)
    print(points)

    L = lagrange_polynomial(points, x)

    print(f"Многочлен Лагранжа: {L}")

    L_diff = sy.diff(L, x)

    f = x ** 2 - sy.cos(sy.pi * x)

    d = take_diff(f, x, n)
    df = take_diff(f, x, m)
    r_1 = d.subs(x, a) - L_diff.subs(x, a)
    r_min = (df.subs(x, a) / sy.factorial(n+1)) * omega(a, b, step, x)

    r_max = (df.subs(x, b) / sy.factorial(n+1)) * omega(a, b, step, x)

    print(f"Производная многочлена Лагранжа: {L_diff}")

    print(L.subs(x, a))

    print(func(a))
    print('----------------------------------------------------------')
    print('Значения дифференциала полинома Лагранжа и дифференциала функции соответственно: ')
    print('----------------------------------------------------------')
    print(L_diff.subs(x, a))
    print(d.subs(x, a))
    print('----------------------------------------------------------')
    print('Значение остаточного члена, верхняя и нижняя оценки: ')
    print('----------------------------------------------------------')
    print(r_1)
    print(r_min)
    print(r_min.subs(x, a))
    print(r_max)
    print(r_max.subs(x, b))



if __name__ == '__main__':
    main()