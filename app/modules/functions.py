import sympy as sy
from sympy.parsing.latex import parse_latex
import inspect
import pickle

def find_symbol(expr: sy.Expr, target: str):
    '''Find a specific symbol in a expression. Function to search for a symbol in a expression'''
    symb_lst = list(expr.free_symbols)
    names_lst = []
    if len(symb_lst) > 0:
        for symbol in expr.free_symbols:

            names_lst.append(symbol.name)

            if symbol.name in target:
                return symbol

            else:
                pass

        return sy.symbols(target, real=True, positive=True)
    # symbol not found, create the same
    else:
        return sy.symbols(target, real=True, positive=True)


# Unify the symbols of different expressions
def unify_symbols(expr1, expr2):
    symbols_expr1 = set([x.name for x in expr1.free_symbols])
    symbols_expr2 = set([x.name for x in expr2.free_symbols])

    common_symbols = symbols_expr1.intersection(symbols_expr2)
    # print(common_symbols)

    symbol_mapping = {}
    for symbol in common_symbols:
        new_symbol = sy.Symbol(symbol)
        symbol_mapping[symbol] = new_symbol

    # print(symbol_mapping)
    expr1 = expr1.subs(symbol_mapping)
    expr2 = expr2.subs(symbol_mapping)

    return expr1, expr2


def save_obj(object, file_name):
    '''Save an pickle object into the disk'''
    with open(file_name, 'wb') as file:
        pickle.dump(object, file)
    print(f"Objeto salvo em {file_name}")


def load_obj(file_name):
    '''Loads a pickle object from the disk'''
    with open(file_name, 'rb') as file:
        object = pickle.load(file)
    return object


# Eq 2.36
def somm_corr():

    expr1 = r"\frac{\pi x}{1-e^{-\pi x}}"
    expr2 = r"\prod_{b = 1}^{l} (1 + \frac{x^2}{4b^2} )"

    fact1 = parse_latex(expr1)
    fact2 = parse_latex(expr2)

    result = fact1*fact2

    return result


def dln(np):
    expr = []
    expr.append(r"\frac{n! (2l+2n+1)!}{(l+n)!}")
    expr.append(
        r"\sum_{j = 0}^{2n} \frac{(-2)^j (l+j)!}{j! (2n-j)! (2l+j+1)!}")
    expr.append(
        r"\left[\prod_{b = l+1}^{l+j}\left(1 + \frac{i x}{2b}\right)\right]")

    fact = []
    result = 1

    for i, exp in enumerate(expr):
        fact.append(parse_latex(exp))
        result *= fact[i]

    nnew = sy.Symbol(f'{np}')
    result = result.replace(find_symbol(result, 'n'), nnew)
    return result


def get_py_src(expression: sy.Expr, func_name: str):
    '''Transform the expressions into python functions automatically'''

    old_symbols = list(expression.free_symbols)
    new_symbols = [str(sym).replace('\\', '_').replace('{', '').replace(
        '}', '') for sym in list(tuple((expression.free_symbols)))]

    new_expression = expression.subs(tuple(zip(old_symbols, new_symbols)))
    py_func = sy.lambdify(
        list(new_expression.free_symbols), new_expression, 'numpy')

    lambdfy_str = inspect.getsource(py_func)
    lambdfy_str = lambdfy_str.replace('_lambdifygenerated', func_name)
    lambdfy_str = lambdfy_str.replace('sqrt', 'np.sqrt')

    return lambdfy_str
