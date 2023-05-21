def process_colormap_file(colormap_file):
    colormap = {}
    with open(colormap_file, 'r') as f:
        for line in f:
            site, color = line.split()
            colormap[site] = color

    return colormap

def get_colormap(sites):
    return {s: (i+1) for i, s in enumerate(sites)}

