import argparse
import pkg_resources

def run():
    version = pkg_resources.get_distribution('c_p_reach').version
    parser = argparse.ArgumentParser(f'c_p_reach {version:s}')
    parser.add_argument("casadi_model")
    args = parser.parse_args()


    try:
        with open(args.casadi_model, 'r') as f:
            script=f.read()
    except FileNotFoundError as e:
        print(e)
        exit(1)

    print(script)

    local_scope = {}
    exec(script, {}, local_scope)
    globals().update(local_scope)
    locals().update(local_scope)
    print('local scope keys', local_scope.keys())
    print('local keys', locals().keys())
    print(Integrator)
    myint = Integrator()
    print(myint)

