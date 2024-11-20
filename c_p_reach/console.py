import argparse
import pkg_resources

def run():
    version = pkg_resources.get_distribution('c_p_reach').version
    parser = argparse.ArgumentParser(f'c-p-reach {version:s}')
    parser.add_argument("casadi_model")
    args = parser.parse_args()

    try:
        with open(args.casadi_model, 'r') as f:
            script=f.read()
    except FileNotFoundError as e:
        print(e)
        exit(1)

    exec(script)
