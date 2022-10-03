import pickaxe_generic as pg

engine = pg.create_engine()
network = engine.network_from_file("network")
print(network)
