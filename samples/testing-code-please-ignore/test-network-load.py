import pickaxe_generic as pg

engine = pg.create_engine()
network_loaded = engine.network_from_file("network")
print(list(network_loaded.mols))
