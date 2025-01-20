"""
A file to configure calling in the prompt line

python -m reedsolomon -e "Hello Word!"

reedsolomon -e "Hello Word!
"""


from pathlib import Path

import sys, os


def main():
    print("For the moment, data can't be encoded.")

    try:
        args: list[str] = [val for i, val in enumerate(sys.argv[1:]) if not i%2]
        vals: list[str] = [arg for i, arg in enumerate(sys.argv[1:]) if i%2]
        dico: dict[str, str] = {arg:val for arg, val in zip(args, vals)}

    except:
        print("Arguments aren't in the right order")
        return 1

    if not '-f' in sys.argv:
        name = "reedsolomon_data.txt"

    else:
        name = dico['-f']

    path = Path(os.path.dirname(sys.argv[0])) / Path(name)

    with path.open('w', encoding='UTF-8') as f:
        data = dico['-e']

        f.write(data)

    print(f"Path : {path}\n{data}")

if __name__ == '__main__':
    sys.exit(main())