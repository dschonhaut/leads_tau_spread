# Write some plotting functions in matplotlib and seaborn.

from matplotlib import pyplot as plt
import seaborn as sns

def plot_scatter(data,
    x_col,
    y_col,
    hue=None,
    palette=None,
    size=None,
    fit_reg=False,
    **kws
    ):
                 
    """Plot a scatter plot with seaborn.

    Options
    -------
    - Can add a line of best fit
    - Can plot groups of points in multiple colors, or with different shapes.

    Parameters
    ----------
    data : pandas.DataFrame
        Data to plot.
    x_col : str
        Column name for x-axis.
    y_col : str
        Column name for y-axis.
    hue : str, optional
        Column name for grouping data.
    palette : str, optional
        Color palette to use for plotting.
    fit_reg : bool, optional
        Plot a line of best fit.
    **kws : optional
        Additional keyword arguments to pass to seaborn.scatterplot.
    """
    fig, ax = plt.subplots()
    ax = np.ravel(ax)

    iax = 0
    _ax = ax[iax]
    sns.scatterplot(
        data=data,
        x=x_col,
        y=y_col,
        hue=hue,
        palette=palette,
        size=size,
        fit_reg=fit_reg,
        ax=_ax,
        **kws
        )
    
    # Format the axis ticks and ticklabels.
    _ax.tick_params(axis='both', which='major', labelsize=12)
    _ax.tick_params(axis='both', which='minor', labelsize=12)
    _ax.xaxis.set_tick_params(rotation=45)
    _ax.yaxis.set_tick_params(rotation=45)

    # Format the axis labels.
    _ax.set_xlabel(x_col, fontsize=14)
    _ax.set_ylabel(y_col, fontsize=14)

    # Format the legend.
    if hue is not None:
        handles, labels = _ax.get_legend_handles_labels()
        _ax.legend(handles=handles[1:], labels=labels[1:])

    fig.tight_layout()

    return fig, ax
    