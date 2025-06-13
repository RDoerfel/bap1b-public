import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path

CM = 2.54

### Classes
class Plot(mpl.figure.Figure):
    def __init__(self, *args, **kwargs) -> mpl.figure.Figure:
        super().__init__(*args, **kwargs)

    def set_style(self, offleft=5, offbottom=5, spinewidth=1.4) -> None:
        """ set style for plots (despine, grid, ticks)
        args:
            offleft: offset left
            offbottom: offset bottom
            spinewidth: width of spines
        """
        sns.despine(self,top=True, right=True, left=False, bottom=False, offset={'left': offleft, 'bottom': offbottom})

        for ax in self.axes:
            ax.grid(axis='y', color='C7', linestyle='--', lw=.8)
            ax.tick_params(which='major', direction='out', length=3, width=spinewidth, bottom=True, left=True)
            ax.tick_params(which='minor', direction='out', length=2, width=spinewidth/2, bottom=True, left=True)
            plt.setp(ax.spines.values(), linewidth=spinewidth)

        self.tight_layout()


    def save(self, filename:str, path:Path=None, **kwargs) -> None:
        """ save figure
        args:
            filename: filename
            path: path to save figure
        """
        if path is None:
            path = Path.cwd()
        self.savefig(path / filename, bbox_inches='tight',dpi=300, **kwargs) 


### functions
def set_axes_description(ax, xlim:tuple=None, ylim:tuple=None, xlable:str=None, ylable:str=None, title:str=None, xticks:list=None, yticks:list=None) -> plt.Axes:
    ax.set_xlabel(xlable)
    ax.set_ylabel(ylable)
    ax.set_title(title)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    if xticks:
        ax.set_xticks(xticks)
    if yticks:
        ax.set_yticks(yticks)
    return ax

def scatter_plot(ax, data, x_col, y_col, scatter_color, line_color, label) -> plt.Axes:
    sns.scatterplot(x=x_col, y=y_col, data=data, ax=ax, color=scatter_color, label=label, legend=False, s=15, alpha=0.3)
    sns.regplot(x=x_col, y=y_col, data=data, ax=ax, scatter=False, color=line_color, line_kws={'linewidth':1.5})
    return ax

def violinplot(ax, data, x_col, y_col, color_dict):
    sns.violinplot(x=x_col, y=y_col, data=data, ax=ax, palette=color_dict, alpha=0.3, linewidth=0.5)
    return ax

def set_rc_params(fontfamily:str=None, small=8, medium=10, big=12) -> None:
    """ set rc parameters for plots 
    args:
        fontfamily: font family
        small: fontsize for small text
        medium: fontsize for medium text
        big: fontsize for big text
    """
    mpl.rcParams.update()
    if fontfamily is not None:
        plt.rc('font', family=fontfamily)
    plt.rc('font', size=medium)         # controls default text sizes
    plt.rc('axes', titlesize=big)       # fontsize of the axes title
    plt.rc('axes', labelsize=medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small)    # legend fontsize
    plt.rc('savefig', dpi=300)          # figure resolution
    plt.rc('pdf', fonttype=42)
    plt.rc('ps', fonttype=42)

def get_figures(rows:int, cols:int, unit:str, figwidth:float, figheight:float, sharex=True,sharey=True) -> list[Plot, plt.Axes]:
    """ Get figure and axes object
    args:
        rows: number of rows
        cols: number of columns
        unit: unit of figure size (cm or inch)
        figwidth: figure width
        figheight: figure height
        sharex: share x axis
    returns: figure object, axes object"""
    if unit == 'cm':
        figsize = (figwidth/CM, figheight/CM)
    elif unit == 'inch':
        figsize = (figwidth, figheight)
    else:
        raise ValueError(f'unit {unit} not supported')
    plot = Plot(figsize=figsize)
    axs = plot.subplots(nrows=rows, ncols=cols, sharex=sharex,sharey=sharey)
    return plot, axs

