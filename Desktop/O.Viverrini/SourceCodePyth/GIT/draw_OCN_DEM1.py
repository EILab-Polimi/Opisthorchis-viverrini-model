import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def draw_OCN_DEM1(OCN, X, 
             Draw_Borders=True, 
             Draw_River=True, 
             Draw_Frame=True,
             Borders_Color='black', 
             River_Color='#2776ea',
             cmap=None,
             alpha=0.6,
             vmin=None,
             vmax=None,
             cbar_label='',
             cbar_label_size=12,
             cbar_tick_size=10,
             show_colorbar=True,
             binary=False):
    """
    Python version of MATLAB draw_OCN function.
    Replicates the same behavior using matplotlib.
    """

    X = np.array(X)

    # If X is per-node, map it onto grid
    if X.shape == (OCN['nNodes'],) or X.shape == (OCN['nNodes'], 1):
        X = assignSC(OCN['MSC'], X)

    if X.shape == OCN['FD']['A'].shape:
        Xtemp = X.copy()
        X = np.zeros_like(OCN['MSC'])
        for i in range(len(Xtemp)):
            xi = OCN['geometry']['XI'][i] - 1
            yi = OCN['geometry']['YI'][i] - 1
            X[xi, yi] = Xtemp.flat[i]

    plt.figure(figsize=(8, 6))

    if not binary:
        if not np.all(np.isnan(X)):

            if cmap is None:
                cmap = plt.cm.pink_r  # default identico a prima

            im = plt.imshow(X.T, cmap=cmap, alpha=alpha, origin='lower')

            # scaling controllabile da fuori
            if vmin is None:
                vmin = np.nanmin(X)
            if vmax is None:
                vmax = np.nanmax(X)

            im.set_clim(vmin, vmax)

            # --- colorbar ---
            if show_colorbar:
                cb = plt.colorbar(im)

                if cbar_label:
                    cb.set_label(cbar_label, fontsize=cbar_label_size)

                cb.ax.tick_params(labelsize=cbar_tick_size)
    else:
        plt.imshow(X.T, origin='lower', alpha=0.4,
                   cmap=(plt.cm.Blues if np.all(X == 0)
                         else ListedColormap(cmap)))

    # --- Draw borders ---
    if Draw_Borders:
        draw_borders(OCN['MSC'], Borders_Color)

    # --- Draw rivers ---
    if Draw_River:
        for i in range(len(OCN['FD']['A'])):
            if i != OCN['outlet'] - 1:
                if OCN['FD']['A'].flat[i] >= OCN['thrA']:
                    x0 = OCN['FD']['X'].flat[i] / OCN['cellsize']
                    y0 = OCN['FD']['Y'].flat[i] / OCN['cellsize']
                    x1 = OCN['FD']['X'].flat[int(OCN['FD']['downNode'][i]-1)] / OCN['cellsize']
                    y1 = OCN['FD']['Y'].flat[int(OCN['FD']['downNode'][i]-1)] / OCN['cellsize']
                    lw = 0.9 * (OCN['FD']['A'].flat[i] * 1e-1 * 1e-8) ** 0.3
                    plt.plot([x0, x1], [y0, y1],
                             color=(0.1451, 0.1176, 0.5412),
                             linewidth=lw,
                             solid_capstyle='round')
                    
    # --- ✅ Draw reservoirs ---
    if 'ReservoirsPixel' in OCN and np.any(np.array(OCN['ReservoirsPixel']) > 0):
        for i in range(len(OCN['ReservoirsPixel'])):
            if OCN['ReservoirsPixel'][i] > 0:
                x = OCN['FD']['X'][i] / OCN['cellsize']
                y = OCN['FD']['Y'][i] / OCN['cellsize']
                plt.plot(x, y, '^', markersize=10,
                         markerfacecolor='r',
                         markeredgecolor='k')

    # --- Draw frame ---
    if Draw_Frame:
        nx, ny = OCN['MSC'].shape
        nx= nx-1
        ny = ny -1
        plt.plot([-0.5, -0.5], [-0.5, ny+0.5], color='k')
        plt.plot([nx+0.5, nx+0.5], [-0.5, ny+0.5], color='k')
        plt.plot([-0.5, nx + .5], [-0.5, -0.5], color='k')
        plt.plot([-0.5, nx + .5], [ny + .5, ny + .5], color='k')

    plt.axis('equal')
    plt.axis('off')
    plt.tight_layout()
    plt.show()


# --- Support function as before ---
def draw_borders(MSC, color='black'):
    nx, ny = MSC.shape
    for i in range(nx-1):
        for j in range(ny-1):
            if MSC[i, j] != MSC[i+1, j]:
                #plt.plot([i + 1, i + 1], [j, j + 1], color=color, linewidth=0.5)
                plt.plot([i + 1-0.5, i + 1-0.5], [j-0.5, j + 1-0.5], color=color, linewidth=0.5)
            if MSC[i, j] != MSC[i, j+1]:
                #plt.plot([i, i + 1], [j + 1, j + 1], color=color, linewidth=0.5)
                plt.plot([i-0.5, i + 1-0.5], [j + 1-0.5, j + 1-0.5], color=color, linewidth=0.5)


def assignSC(MSC, values):
    out = np.zeros_like(MSC)
    for i in range(1, int(np.max(MSC)) + 1):
        out[MSC == i] = values[i - 1]
    return out
