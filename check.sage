
from itertools import count
from itertools import combinations
from itertools import chain
from random import sample

def to_linear_program(polyhedron, solver=None):

    from sage.rings.rational_field import QQ
    R = polyhedron.base_ring()
    if (solver is not None and
        solver.lower() == 'ppl' and
        R.is_exact() and (not R == QQ)):
        raise NotImplementedError('Cannot use PPL on exact irrational data.')

    from sage.numerical.mip import MixedIntegerLinearProgram
    p = MixedIntegerLinearProgram(solver=solver)
    x = p.new_variable(integer=True, nonnegative=False)

    for ineqn in polyhedron.inequalities_list():
        b = -ineqn.pop(0)
        p.add_constraint(p.sum([x[i]*ineqn[i] for i in range(len(ineqn))]) >= b)

    for eqn in polyhedron.equations_list():
        b = -eqn.pop(0)
        p.add_constraint(p.sum([x[i]*eqn[i] for i in range(len(eqn))]) == -b)

    return p


def powerset(s):
    return map( frozenset , chain.from_iterable(combinations(s, r) for r in range(1,len(s)+1)) )

def partitions(s):
    for A in powerset(s) :
        yield A , s - A

def check ( k ) :

    for v in count(k) :

        for m in range( 1 , v + 1 ) :

            a = v // m
            if a != v / m:
                continue

            for n in range( 1 , a + 1 ) :

                l = a // n
                ncells = m * n * l

                if ncells != v or m > n or n > l or l != a / n or m == n == 1:
                    continue

                print( v , m , n , l )

                def xyz ( i ) :

                    z = i % l
                    i //= l
                    y = i % n
                    i //= n
                    x = i

                    return x , y , z

                for C in combinations( range(ncells) , k ) :

                    for A , B in partitions( frozenset( C ) ) :

                        P = Polyhedron( map( xyz , A ) )
                        Q = Polyhedron( map( xyz , B ) )

                        I = P.intersection(Q)

                        p = to_linear_program(I,solver="GLPK")

                        try:

                            p.solve()
                            break

                        except:
                            pass

                    else:
                        yield C






if __name__ == '__main__' :

    import sys
    for weird_animal in check( int(sys.argv[1]) ) :

        print('BINGO' , weird_animal)
