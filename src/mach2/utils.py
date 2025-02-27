def process_colormap_file(colormap_file):
    """Takes as input a colormap file and returns a colormap (dictionary)."""
    colormap = {}
    with open(colormap_file, 'r') as f:
        for line in f:
            site, color = line.split()
            colormap[site] = color
    return colormap

def get_colormap(locations):
    """Takes as input a collection of locations and returns a colormap (dictionary)."""
    return {s: (i+1) for i, s in enumerate(locations)}

