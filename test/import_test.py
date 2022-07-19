from pkgutil import iter_modules
from importlib import import_module
import pickaxe_generic

pkg_name = pickaxe_generic.__name__
pkg_path = pickaxe_generic.__path__

for _, modname, _ in iter_modules(pkg_path,f"{pkg_name}."):
    module = import_module(modname, "pickaxe_generic")
    print(f"Imported {module}")