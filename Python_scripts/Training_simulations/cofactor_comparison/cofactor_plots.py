import pandas as pd
import matplotlib.pyplot as plt
import pickle 

model_names = ["IFO0880", "IFO0880_jsb", 'iRhtoC', 'rhtoGEM']
reactions = set()
for m in model_names:
    with open (f"Python_scripts\\cofactor_comparison\\{m}_cofactors.pk", "rb") as f:
        df = pickle.load(f)
        reactions.update(df[0].index)
        
print(reactions)        

# model_names = ("Adelie", "Chinstrap", "Gentoo")
# penguin_means = {
#     'Bill Depth': (18.35, 18.43, 14.98),
#     'Bill Length': (38.79, 48.83, 47.50),
#     'Flipper Length': (189.95, 195.82, 217.19),
# }

# x = range(len(model_names))  # the label locations
# width = 1/(len(model_names)+1)  # the width of the bars
# multiplier = 0

# fig, ax = plt.subplots(layout='constrained')

# for i in range(len(model_names)):
#     offset = width * multiplier
#     rects = ax.bar(x + offset, measurement, width, label=attribute)
#     ax.bar_label(rects, padding=3)
#     multiplier += 1