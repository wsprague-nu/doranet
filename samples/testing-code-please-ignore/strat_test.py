import pickaxe_generic as pg

engine = pg.create_engine()
network = engine.new_network()
strat = engine.strat.cartesian(network)
