from operator import add

from pickaxe_generic.utils import logreduce

for i in range(1, 10):
    print(f"{logreduce(add,range(i))}")
