



import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import flare_lf.evo as evo



model = 'flares_schechter'

m = evo.read(model, scheme = 'interp')

print(m)
print(m.t.meta['ads'])


print(m.p(5.0))



model = 'flares_binned'

m = evo.read(model)

print(m.p(5.0))
