import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
import shutil
from matplotlib.colors import ListedColormap, Normalize
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from shapely.geometry import Polygon
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

DATA_DIR = '../../data/'
RESULTS_DIR = '../../results/'
FIGURES_DIR = '../../figures/'

METHODS_ORDER = ['CPDA', 'Terminal-Dominated Matching', 'Multi-Conn Heuristic', 'Coordinated BH']
METHOD_MARKERS = {
    'CPDA': 'o',
    'Terminal-Dominated Matching': 'D',
    'Multi-Conn Heuristic': 's',
    'Coordinated BH': '^'
}
METHOD_COLORS = {
    'CPDA': '#1f77b4',
    'Terminal-Dominated Matching': '#2ca02c',
    'Multi-Conn Heuristic': '#ff7f0e',
    'Coordinated BH': '#9467bd',
    'CPDA (Without Fragmentation)': '#d62728',
}

def save_pdf_svg(path):
    plt.savefig(path, format='pdf', bbox_inches='tight')
    plt.savefig(path.replace('.pdf', '.svg'), format='svg', bbox_inches='tight')

def method_columns_to_long(df, prefix, value_name):
    rename_map = {
        f'{prefix}_CPDA': 'CPDA',
        f'{prefix}_TDM': 'Terminal-Dominated Matching',
        f'{prefix}_MultiConn': 'Multi-Conn Heuristic',
        f'{prefix}_CBH': 'Coordinated BH'
    }
    long_df = df.rename(columns=rename_map).melt(
        id_vars=['cell_id'],
        value_vars=[m for m in METHODS_ORDER if m in df.rename(columns=rename_map).columns],
        var_name='method',
        value_name=value_name
    )
    long_df[value_name] = pd.to_numeric(long_df[value_name], errors='coerce')
    return long_df.dropna(subset=[value_name])

def style_axis(ax):
    ax.grid(True, axis='y', linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color('k')

def draw_boxplot(ax, arrays, positions, widths, colors, showmeans=False):
    bp = ax.boxplot(
        arrays,
        positions=positions,
        widths=widths,
        patch_artist=True,
        showfliers=False,
        showmeans=showmeans,
        medianprops={'color': 'k', 'linewidth': 1.0},
        whiskerprops={'color': 'k', 'linewidth': 0.8},
        capprops={'color': 'k', 'linewidth': 0.8},
        meanprops={'marker': 'o', 'markerfacecolor': 'none',
                   'markeredgecolor': 'k', 'markersize': 4}
    )
    for box, color in zip(bp['boxes'], colors):
        box.set(facecolor=color, edgecolor='k', alpha=0.68, linewidth=0.8)
    return bp

def plot_grouped_box_by_method(df, group_col, group_order, group_labels, out_pdf,
                               ylabel, xlabel, ylim=None, legend_loc='best'):
    df = df.dropna(subset=[group_col, 'method', 'value']).copy()
    df['value'] = pd.to_numeric(df['value'], errors='coerce')
    df = df.dropna(subset=['value'])

    present_methods = set(df['method'].unique())
    methods = [m for m in METHODS_ORDER if m in present_methods]
    groups = [g for g in group_order if g in set(df[group_col].unique())]
    if not methods or not groups:
        raise ValueError(f'No plottable data for {out_pdf}')

    x = np.arange(len(groups))
    width = min(0.18, 0.72 / len(methods))
    offsets = (np.arange(len(methods)) - (len(methods) - 1) / 2.0) * width

    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    legend_handles = []
    for j, method in enumerate(methods):
        color = METHOD_COLORS.get(method, plt.cm.tab10(j))
        arrays = [
            df[(df[group_col] == g) & (df['method'] == method)]['value'].values
            for g in groups
        ]
        xpos = x + offsets[j]
        draw_boxplot(ax, arrays, xpos, width * 0.62, [color] * len(arrays))
        legend_handles.append(Patch(facecolor=color, edgecolor='k', alpha=0.68, label=method))

    ax.set_xticks(x)
    ax.set_xticklabels([group_labels.get(g, str(g)) for g in groups], fontsize=18)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    if ylim is not None:
        ax.set_ylim(ylim)
    style_axis(ax)
    ax.legend(handles=legend_handles, fontsize=14, loc=legend_loc, ncol=1, framealpha=0.88)
    fig.tight_layout()
    save_pdf_svg(out_pdf)
    plt.close(fig)

def plot_cell_metrics_box_panels(metric_frames, out_pdf):
    methods = [m for m in METHODS_ORDER if any(m in set(df['method'].unique()) for df, _, _ in metric_frames)]
    colors = [METHOD_COLORS.get(m, plt.cm.tab10(i)) for i, m in enumerate(methods)]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.8), sharex=True)

    for ax, (df, ylabel, title) in zip(axes, metric_frames):
        arrays = [df[df['method'] == m][ylabel].values for m in methods]
        positions = np.arange(1, len(methods) + 1)
        draw_boxplot(ax, arrays, positions, 0.46, colors, showmeans=False)
        ax.set_title(title, fontsize=18)
        ax.set_ylabel(ylabel, fontsize=18)
        ax.set_xticks([])
        if 'Unused' in ylabel:
            ax.set_yscale('symlog', linthresh=0.1)
            ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
        style_axis(ax)

    handles = [
        Patch(facecolor=METHOD_COLORS.get(m, plt.cm.tab10(i)), edgecolor='k', alpha=0.68, label=m)
        for i, m in enumerate(methods)
    ]
    fig.legend(handles=handles, loc='upper center', ncol=4, fontsize=14,
               framealpha=0.88, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    save_pdf_svg(out_pdf)
    plt.close(fig)

def plot_cell_metrics_grouped_bars(df_rate, df_unmet, df_unused, out_pdf):
    cell_ids = pd.to_numeric(df_rate['cell_id'], errors='coerce').astype(int).values
    x = np.arange(len(cell_ids))

    method_columns = [
        ('CPDA', 'CPDA'),
        ('Terminal-Dominated Matching', 'TDM'),
        ('Multi-Conn Heuristic', 'MultiConn'),
        ('Coordinated BH', 'CBH')
    ]
    demand_color = '#7f7f7f'

    panels = [
        (
            'Rate',
            'Average Rate (Mbps)',
            [(df_rate['D_bar'].values, 'Mean Demand', demand_color)] +
            [(df_rate[f'Rbar_{suffix}'].values, label, METHOD_COLORS[label])
             for label, suffix in method_columns],
            False
        ),
        (
            'Unmet',
            'Unmet Demand (Mbps)',
            [(df_unmet[f'Unmetbar_{suffix}'].values, label, METHOD_COLORS[label])
             for label, suffix in method_columns],
            False
        ),
        (
            'Unused',
            'Unused Rate (Mbps)',
            [(df_unused[f'Unusedbar_{suffix}'].values, label, METHOD_COLORS[label])
             for label, suffix in method_columns],
            True
        )
    ]

    fig, axes = plt.subplots(3, 1, figsize=(17.2, 10.2), sharex=True)
    legend_handles = {}
    for ax, (title, ylabel, series, use_symlog) in zip(axes, panels):
        n_bars = len(series)
        group_width = 0.82
        bar_width = group_width / n_bars
        offsets = (np.arange(n_bars) - (n_bars - 1) / 2.0) * bar_width

        for offset, (values, label, color) in zip(offsets, series):
            values = pd.to_numeric(pd.Series(values), errors='coerce').fillna(0).values
            ax.bar(
                x + offset,
                values,
                width=bar_width * 0.92,
                color=color,
                alpha=0.78,
                edgecolor='k',
                linewidth=0.45,
                label=label
            )
            legend_handles.setdefault(label, Patch(facecolor=color, edgecolor='k', alpha=0.78, label=label))

        ax.set_title(title, fontsize=18)
        ax.set_ylabel(ylabel, fontsize=18)
        ax.set_xlim(-0.7, len(cell_ids) - 0.3)
        ax.set_ylim(bottom=0)
        if use_symlog:
            ax.set_yscale('symlog', linthresh=0.1)
            ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
        style_axis(ax)

    axes[-1].set_xlabel('Cell ID', fontsize=18)
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(cell_ids, fontsize=11)

    ordered_labels = ['Mean Demand'] + [label for label, _ in method_columns]
    handles = [legend_handles[label] for label in ordered_labels if label in legend_handles]
    fig.legend(handles=handles, loc='upper center', ncol=5, fontsize=14,
               framealpha=0.88, bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    save_pdf_svg(out_pdf)
    plt.close(fig)

def plot_method_summary_bars(long_df, value_col, out_pdf, ylabel, demand_series=None):
    long_df = long_df.dropna(subset=['method', value_col]).copy()
    methods = [m for m in METHODS_ORDER if m in set(long_df['method'].unique())]
    means = [long_df[long_df['method'] == m][value_col].mean() for m in methods]
    stds = [long_df[long_df['method'] == m][value_col].std(ddof=1) for m in methods]
    colors = [METHOD_COLORS.get(m, plt.cm.tab10(i)) for i, m in enumerate(methods)]

    fig, ax = plt.subplots(figsize=(9, 6))
    x = np.arange(len(methods))
    ax.bar(x, means, yerr=stds, capsize=5, color=colors, alpha=0.78,
           edgecolor='k', linewidth=0.8)
    if demand_series is not None:
        ax.axhline(pd.to_numeric(demand_series, errors='coerce').mean(),
                   linestyle='--', color='k', linewidth=1.2, label='Mean Demand')
        ax.legend(fontsize=14, loc='upper right')
    ax.set_xticks(x)
    labels = [
        'CPDA\n(Without Fragmentation)' if m == 'CPDA (Without Fragmentation)' else m
        for m in methods
    ]
    ax.set_xticklabels(labels, rotation=0 if value_col == 'ASE' else 18,
                       ha='center' if value_col == 'ASE' else 'right',
                       fontsize=16)
    ax.set_ylabel(ylabel, fontsize=18)
    if value_col == 'ASE':
        formatter = mticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-4, -4))
        ax.yaxis.set_major_formatter(formatter)
        ax.yaxis.get_offset_text().set_fontsize(16)
    ax.grid(True, axis='y', linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    fig.tight_layout()
    save_pdf_svg(out_pdf)
    plt.close(fig)

def plot_box_by_method(df, value_col, out_pdf, ylabel, methods_order, colors, linestyles=None):
    df = df.dropna(subset=['method', value_col]).copy()
    df[value_col] = pd.to_numeric(df[value_col], errors='coerce')
    df = df.dropna(subset=[value_col])
    methods = [m for m in methods_order if m in set(df['method'].unique())]
    arrays = [df[df['method'] == m][value_col].values for m in methods]
    fig, ax = plt.subplots(figsize=(8, 5.8))
    box_colors = [colors.get(m, plt.cm.tab10(i)) for i, m in enumerate(methods)]
    draw_boxplot(ax, arrays, np.arange(1, len(methods) + 1),
                 0.26 * 0.68, box_colors, showmeans=False)
    if value_col == 'ASE':
        ax.set_xticks([(len(methods) + 1) / 2])
        ax.set_xticklabels(['30 MHz'], fontsize=18)
        ax.set_xlabel('System Bandwidth', fontsize=18)
        handles = [
            Patch(facecolor=colors.get(m, plt.cm.tab10(i)), edgecolor='k', alpha=0.68,
                  label=('CPDA (Without Fragmentation)' if m == 'CPDA (Without Fragmentation)' else m))
            for i, m in enumerate(methods)
        ]
        ax.legend(handles=handles, loc='lower left', fontsize=14, framealpha=0.88)
    else:
        ax.set_xticklabels(methods, rotation=18, ha='right', fontsize=16)
    ax.set_ylabel(ylabel, fontsize=18)
    if value_col == 'ASE':
        formatter = mticker.ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((-4, -4))
        ax.yaxis.set_major_formatter(formatter)
        ax.yaxis.get_offset_text().set_fontsize(16)
    style_axis(ax)
    fig.tight_layout()
    save_pdf_svg(out_pdf)
    plt.close(fig)

def plot_bandwidth_ase_boxes(df, out_pdf):
    df = df.copy()
    df['ASE_frag'] = pd.to_numeric(df['ASE_frag'], errors='coerce')
    df['ASE_nofrag'] = pd.to_numeric(df['ASE_nofrag'], errors='coerce')
    df = df.dropna(subset=['B_MHz', 'ASE_frag', 'ASE_nofrag'])
    long_df = df.melt(id_vars=['B_MHz'], value_vars=['ASE_frag', 'ASE_nofrag'],
                      var_name='mode', value_name='ASE')
    mode_labels = {'ASE_frag': 'CPDA', 'ASE_nofrag': 'CPDA (Without Fragmentation)'}
    modes = ['ASE_frag', 'ASE_nofrag']
    bandwidths = sorted(long_df['B_MHz'].unique())
    x = np.arange(len(bandwidths))
    width = 0.26
    colors = {
        'ASE_frag': METHOD_COLORS['CPDA'],
        'ASE_nofrag': METHOD_COLORS['CPDA (Without Fragmentation)']
    }

    fig, ax = plt.subplots(figsize=(8, 6))
    handles = []
    for j, mode in enumerate(modes):
        arrays = [long_df[(long_df['B_MHz'] == b) & (long_df['mode'] == mode)]['ASE'].values
                  for b in bandwidths]
        xpos = x + (j - 0.5) * width
        draw_boxplot(ax, arrays, xpos, width * 0.68, [colors[mode]] * len(arrays))
        means = [np.mean(a) for a in arrays]
        ax.plot(xpos, means, linestyle='--' if mode == 'ASE_nofrag' else '-',
                marker='o', markersize=5, markerfacecolor='none',
                color=colors[mode], linewidth=1)
        handles.append(Line2D([0], [0], color=colors[mode],
                              linestyle='--' if mode == 'ASE_nofrag' else '-',
                              marker='o', markerfacecolor='none',
                              linewidth=1, label=mode_labels[mode]))
    ax.set_xticks(x)
    ax.set_xticklabels([f'{int(b)} MHz' for b in bandwidths], fontsize=18)
    ax.set_xlabel('System Bandwidth', fontsize=18)
    ax.set_ylabel('Area Spectral Efficiency (bit/s/Hz/km$^2$)', fontsize=18)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-4, -4), useMathText=True)
    ax.yaxis.get_offset_text().set_fontsize(16)
    style_axis(ax)
    ax.legend(handles=handles, loc='lower right', fontsize=14, framealpha=0.88)
    fig.tight_layout()
    save_pdf_svg(out_pdf)
    plt.close(fig)

def plot_runtime_boxplot(df, out_pdf, method_order=None):
    if method_order is None:
        method_order = [
            'CPDA',
            'Multi-Conn Heuristic',
            'Coordinated BH',
            'Terminal-Dominated Matching'
        ]
    runtime_colors = {
        'Layout Generation': '#7f7f7f',
        'CPDA': METHOD_COLORS['CPDA'],
        'Multi-Conn Heuristic': METHOD_COLORS['Multi-Conn Heuristic'],
        'Coordinated BH': METHOD_COLORS['Coordinated BH'],
        'Terminal-Dominated Matching': METHOD_COLORS['Terminal-Dominated Matching']
    }
    display_labels = {
        'Layout Generation': 'Layout Generation',
        'CPDA': 'Proposed',
        'Multi-Conn Heuristic': 'Multi-Conn Heuristic',
        'Coordinated BH': 'Coordinated BH',
        'Terminal-Dominated Matching': 'Terminal-Dominated Matching'
    }

    df = df.dropna(subset=['method', 'runtime_s']).copy()
    df['runtime_s'] = pd.to_numeric(df['runtime_s'], errors='coerce')
    df = df.dropna(subset=['runtime_s'])
    df = df[df['runtime_s'] > 0]
    methods = [method for method in method_order if method in set(df['method'].unique())]
    if not methods:
        raise ValueError(f'No runtime data are available for {out_pdf}.')

    arrays = [df[df['method'] == method]['runtime_s'].values for method in methods]
    positions = np.arange(1, len(methods) + 1)
    fig_height = 3.8 if len(methods) <= 2 else 5.0
    fig, ax = plt.subplots(figsize=(9.2, fig_height))
    bp = ax.boxplot(
        arrays,
        positions=positions,
        widths=0.18,
        vert=False,
        patch_artist=True,
        showfliers=False,
        medianprops={'linewidth': 0.7},
        whiskerprops={'linewidth': 0.6},
        capprops={'linewidth': 0.6}
    )
    for idx, (box, method) in enumerate(zip(bp['boxes'], methods)):
        color = runtime_colors[method]
        box.set(
            facecolor=mpl.colors.to_rgba(color, 0.24),
            edgecolor='k',
            alpha=1.0,
            linewidth=0.65
        )
        bp['medians'][idx].set_color('k')
        for whisker in bp['whiskers'][2 * idx:2 * idx + 2]:
            whisker.set_color('k')
        for cap in bp['caps'][2 * idx:2 * idx + 2]:
            cap.set_color('k')

    ax.set_yticks(positions)
    ax.set_yticklabels([''] * len(methods))
    ax.tick_params(axis='y', left=False)
    ax.invert_yaxis()
    ax.set_xscale('log')
    ax.xaxis.set_major_locator(mticker.LogLocator(base=10.0))
    ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
    ax.set_xlabel('Runtime (s)', fontsize=18)
    ax.grid(True, axis='x', linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color('k')
    legend_handles = [
        Patch(facecolor=mpl.colors.to_rgba(runtime_colors[method], 0.45),
              edgecolor='k', linewidth=0.65, label=display_labels[method])
        for method in methods
    ]
    fig.legend(handles=legend_handles, loc='upper center', ncol=len(methods),
               fontsize=13, framealpha=0.88, bbox_to_anchor=(0.5, 1.01))
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    save_pdf_svg(out_pdf)
    plt.close(fig)

mpl.rcParams['font.family'] = 'Times New Roman'

# Figure 1: cell demand density and satellite positions.
cell_data = pd.read_csv(DATA_DIR + 'cell_clusters/output_cell_cluster1.csv')
center_data = cell_data[cell_data['region'] == 'center'].reset_index(drop=True)


visibility_data = pd.read_csv(DATA_DIR + 'snapshots_Phase2/visibility_relations_cluster1_config2.csv')
first_snapshot = visibility_data['time'].min()
satellite_positions = visibility_data[visibility_data['time'] == first_snapshot][
    ['satellite_id', 'sat_lat_deg', 'sat_lon_deg', 'sat_height_km']
].drop_duplicates().rename(columns={
    'satellite_id': 'SatelliteIndex',
    'sat_lat_deg': 'Latitude',
    'sat_lon_deg': 'Longitude',
    'sat_height_km': 'Altitude'
})


density_data = pd.read_csv(DATA_DIR + 'relative_density.csv')


relative_density = density_data['relative_density'].values

num_cells = center_data.shape[0]

all_cell_lats = []
all_cell_lons = []
for field in ['vertex1_lat', 'vertex2_lat', 'vertex3_lat', 'vertex4_lat', 'vertex5_lat', 'vertex6_lat']:
    all_cell_lats.extend(center_data[field].to_numpy())
for field in ['vertex1_lon', 'vertex2_lon', 'vertex3_lon', 'vertex4_lon', 'vertex5_lon', 'vertex6_lon']:
    all_cell_lons.extend(center_data[field].to_numpy())

sat_lats = satellite_positions['Latitude'].values
sat_lons = satellite_positions['Longitude'].values

all_lats = np.concatenate((np.array(all_cell_lats), sat_lats))
all_lons = np.concatenate((np.array(all_cell_lons), sat_lons))
minLat, maxLat = np.min(all_lats), np.max(all_lats)
minLon, maxLon = np.min(all_lons), np.max(all_lons)
latMargin = 0.05 * (maxLat - minLat)
lonMargin = 0.05 * (maxLon - minLon)
extent = [minLon - lonMargin, maxLon + lonMargin, minLat - latMargin, maxLat + latMargin]
center_lat = (minLat + maxLat) / 2
center_lon = (minLon + maxLon) / 2

cmap_array = np.column_stack((np.linspace(1, 1, 256),
                              np.linspace(0.8, 0, 256),
                              np.linspace(0.8, 0, 256)))
custom_cmap = ListedColormap(cmap_array)
norm = Normalize(vmin=np.min(relative_density), vmax=np.max(relative_density))

fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection=ccrs.Mercator())
ax.set_extent(extent, crs=ccrs.PlateCarree())

ax.add_feature(cfeature.LAND, facecolor='lightgray')
ax.add_feature(cfeature.OCEAN, facecolor='white')
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')

for i in range(num_cells):
    vertices_lat = [
        center_data.loc[i, 'vertex1_lat'],
        center_data.loc[i, 'vertex2_lat'],
        center_data.loc[i, 'vertex3_lat'],
        center_data.loc[i, 'vertex4_lat'],
        center_data.loc[i, 'vertex5_lat'],
        center_data.loc[i, 'vertex6_lat'],
        center_data.loc[i, 'vertex1_lat']
    ]
    vertices_lon = [
        center_data.loc[i, 'vertex1_lon'],
        center_data.loc[i, 'vertex2_lon'],
        center_data.loc[i, 'vertex3_lon'],
        center_data.loc[i, 'vertex4_lon'],
        center_data.loc[i, 'vertex5_lon'],
        center_data.loc[i, 'vertex6_lon'],
        center_data.loc[i, 'vertex1_lon']
    ]
    poly = Polygon(list(zip(vertices_lon, vertices_lat)))
    fill_color = custom_cmap(norm(relative_density[i]))
    ax.add_geometries([poly], crs=ccrs.PlateCarree(), facecolor=fill_color, edgecolor='black', alpha=0.7)

for k in range(len(sat_lats)):
    ax.plot(sat_lons[k], sat_lats[k], marker='o', markersize=6,
            markerfacecolor='#0050A0', markeredgecolor='black', markeredgewidth = 0.8,
            transform=ccrs.PlateCarree())

sm = plt.cm.ScalarMappable(cmap=custom_cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, orientation='vertical', pad=0.05)
cbar.set_label('Traffic Density')


from matplotlib.lines import Line2D
from matplotlib.patches import Patch

legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Satellite',
           markerfacecolor='#0050A0', markersize=6, markeredgecolor='black', markeredgewidth=0.8),
    Line2D([0], [0], marker='h', color='w', label='Cell',
           markerfacecolor='#FF6666', markersize=15, markeredgecolor='black')
]

plt.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, -0.05),
           fontsize=12, ncol=2, frameon=True)



from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
gl.top_labels = False
gl.right_labels = False

gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

gl.xlabel_style = {'fontname': 'Times New Roman', 'size': 10}
gl.ylabel_style = {'fontname': 'Times New Roman', 'size': 10}

gl.xlocator = mticker.FixedLocator(np.arange(np.floor(minLon/10)*10, np.ceil(maxLon/10)*10+10, 10))
gl.ylocator = mticker.FixedLocator(np.arange(np.floor(minLat/10)*10, np.ceil(maxLat/10)*10+10, 10))





for artist in ax.collections:
    artist.set_rasterized(True)

for artist in ax.images:
    artist.set_rasterized(True)

fig.set_size_inches(6, 5)

# Add a 50 km scale bar.
from pyproj import Geod
_scale_geod = Geod(ellps='WGS84')
_scale_x, _scale_y = ax.transData.inverted().transform(
    ax.transAxes.transform((0.055, 0.065)))
_scale_lon, _scale_lat = ccrs.PlateCarree().transform_point(
    _scale_x, _scale_y, ax.projection)
_scale_low, _scale_high = 0.0, 2.0
for _ in range(50):
    _scale_delta = (_scale_low + _scale_high) / 2
    _scale_distance = _scale_geod.inv(
        _scale_lon, _scale_lat, _scale_lon + _scale_delta, _scale_lat)[2]
    if _scale_distance < 50000:
        _scale_low = _scale_delta
    else:
        _scale_high = _scale_delta
_scale_end_x, _ = ax.projection.transform_point(
    _scale_lon + _scale_delta, _scale_lat, ccrs.PlateCarree())
_scale_tick = 0.012 * (ax.get_ylim()[1] - ax.get_ylim()[0])
ax.plot([_scale_x, _scale_end_x], [_scale_y, _scale_y],
        color='black', linewidth=1.2, zorder=20)
for _scale_tick_x in (_scale_x, _scale_end_x):
    ax.plot([_scale_tick_x, _scale_tick_x],
            [_scale_y - _scale_tick, _scale_y + _scale_tick],
            color='black', linewidth=1.2, zorder=20)
ax.text((_scale_x + _scale_end_x) / 2, _scale_y + 1.6 * _scale_tick,
        '50 km', ha='center', va='bottom', fontsize=10,
        fontname='Times New Roman', color='black', zorder=20)

output_path = FIGURES_DIR + 'figure1.pdf'
plt.savefig(output_path, format='pdf', dpi=300, bbox_inches='tight')




# Figure 2: center and wrap-around cell layout.
fig2_data = cell_data


lat_columns = ['vertex1_lat', 'vertex2_lat', 'vertex3_lat', 'vertex4_lat', 'vertex5_lat', 'vertex6_lat']
lon_columns = ['vertex1_lon', 'vertex2_lon', 'vertex3_lon', 'vertex4_lon', 'vertex5_lon', 'vertex6_lon']

all_lats = np.concatenate([fig2_data[col].to_numpy() for col in lat_columns])
all_lons = np.concatenate([fig2_data[col].to_numpy() for col in lon_columns])
min_lat = np.min(all_lats)
max_lat = np.max(all_lats)
min_lon = np.min(all_lons)
max_lon = np.max(all_lons)

fig, ax = plt.subplots(figsize=(8, 6))

for idx, row in fig2_data.iterrows():
    vertices_lat = [row['vertex1_lat'], row['vertex2_lat'], row['vertex3_lat'],
                    row['vertex4_lat'], row['vertex5_lat'], row['vertex6_lat'], row['vertex1_lat']]
    vertices_lon = [row['vertex1_lon'], row['vertex2_lon'], row['vertex3_lon'],
                    row['vertex4_lon'], row['vertex5_lon'], row['vertex6_lon'], row['vertex1_lon']]

    if row['region'].strip().lower() == 'center':
        fillColor = (0.85, 0.85, 0.85)
        edgeColor = (0.4, 0.4, 0.4)
    else:
        fillColor = (0.95, 0.95, 0.95)
        edgeColor = (0.7, 0.7, 0.7)

    ax.fill(vertices_lon, vertices_lat, color=fillColor, edgecolor=edgeColor, alpha=0.7, linewidth=1)

    if row['region'].strip().lower() == 'center':
        ax.text(row['center_lon'], row['center_lat'], str(int(row['index'])),
                ha='center', va='center', fontsize=10, color=(0.4, 0.4, 0.4), fontweight='bold',
                fontname='Times New Roman')
    else:
        if int(row['index']) in [1, 22, 25, 28, 31, 34, 37]:
            ax.text(row['center_lon'], row['center_lat'], str(int(row['index'])),
                    ha='center', va='center', fontsize=10, color='k', fontweight='bold',
                    fontname='Times New Roman')


ax.set_xticks([])
ax.set_yticks([])


output_path_fig2 = FIGURES_DIR + 'figure2.pdf'
plt.savefig(output_path_fig2, format='pdf', bbox_inches='tight')


mpl.rcParams['font.family'] = 'Times New Roman'

# Figures 3-6: beam illumination patterns and marginal occupancy profiles.
connCount_C = pd.read_csv(RESULTS_DIR + 'figure3_connCount_C.csv', header=None).values
fig3_info = pd.read_csv(RESULTS_DIR + 'figure3_info.csv')
N_period_fig3 = int(fig3_info['N_period'].iloc[0])
num_cells_fig3 = int(fig3_info['num_cells'].iloc[0])

connCount_C_ultimate = pd.read_csv(RESULTS_DIR + 'figure4_connCount_C_ultimate.csv', header=None).values
fig4_info = pd.read_csv(RESULTS_DIR + 'figure4_info.csv')
N_period_fig4 = int(fig4_info['N_period'].iloc[0])
num_cells_fig4 = int(fig4_info['num_cells'].iloc[0])

connCount_C_ultimate_MC = pd.read_csv(
    RESULTS_DIR + 'figure5_connCount_C_ultimate_MC.csv',
    header=None
).values
fig5_info = pd.read_csv(RESULTS_DIR + 'figure5_info.csv')
N_period_fig5 = int(fig5_info['N_period'].iloc[0])
num_cells_fig5 = int(fig5_info['num_cells'].iloc[0])

connCount_C_ultimate_TDM = pd.read_csv(
    RESULTS_DIR + 'figure6_connCount_C_ultimate_TDM.csv',
    header=None
).values
fig6_info = pd.read_csv(RESULTS_DIR + 'figure6_info.csv')
N_period_fig6 = int(fig6_info['N_period'].iloc[0])
num_cells_fig6 = int(fig6_info['num_cells'].iloc[0])


_BEAM_CONNECTION_COLORS = [
    '#3787C0', '#52B8CD', '#20B5AA', '#70BC70',
    '#B6CB45', '#F2CE37', '#F59A27', '#E52B2F'
]
_BEAM_CONNECTION_CMAP = ListedColormap(['#FFFFFF'] + _BEAM_CONNECTION_COLORS)
_BEAM_CONNECTION_CMAP.set_under('white')
_BEAM_CONNECTION_NORM = mpl.colors.BoundaryNorm(
    np.arange(-0.5, 9.5, 1), _BEAM_CONNECTION_CMAP.N)


def _beam_constant_positive_runs(row):
    """Yield start slot, duration, and count for constant nonzero runs."""
    start = None
    value = 0
    for idx, current in enumerate(np.r_[row, 0]):
        if start is None and current > 0:
            start = idx
            value = int(current)
        elif start is not None and current != value:
            yield start, idx - start, value
            start = idx if current > 0 else None
            value = int(current)


def _style_beam_pattern_axis(ax, n_slots, n_cells):
    """Apply the paper's axis style without a background matrix grid."""
    ax.set_xlim(0.5, n_slots + 0.5)
    ax.set_ylim(0.5, n_cells + 0.5)
    ax.set_xticks(np.arange(1, n_slots + 1, 10))
    ax.set_yticks(np.arange(1, n_cells + 1, 5))
    ax.set_xlabel('Time Slot', fontname='Times New Roman', fontsize=18, labelpad=3)
    ax.set_ylabel('Global Cell Index', fontname='Times New Roman', fontsize=18)
    ax.tick_params(axis='both', labelsize=12)
    ax.grid(False)
    ax.set_facecolor('white')
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(0.8)


def _draw_beam_gantt(ax, matrix, method_color):
    """Draw contiguous illumination ribbons with count-dependent intensity."""
    shades = _BEAM_CONNECTION_COLORS
    for cell, row in enumerate(matrix, start=1):
        for start, duration, overlap_count in _beam_constant_positive_runs(row):
            ax.broken_barh(
                [(start + 0.55, max(duration - 0.10, 0.18))],
                (cell - 0.34, 0.68),
                facecolors=shades[overlap_count - 1],
                edgecolors='none',
                linewidth=0
            )


def _draw_active_cell_profile(ax, matrix, method_color, n_slots, active_y_max):
    """Plot the number of active cells in every scheduling slot."""
    active_cells = np.count_nonzero(matrix, axis=0)
    slots = np.arange(1, n_slots + 1)
    for reference_level in (10, 20):
        ax.axhline(reference_level, color='#D9D9D9', linewidth=0.55, zorder=0)
    ax.plot(slots, active_cells, color='black', linewidth=1.25)
    ax.set_xlim(0.5, n_slots + 0.5)
    ax.set_ylim(0, active_y_max)
    ax.set_yticks(np.arange(0, active_y_max + 1, 10))
    ax.set_ylabel('Active\nCells', fontsize=18, fontname='Times New Roman',
                  rotation=0, ha='right', va='center', labelpad=10)
    ax.tick_params(axis='x', bottom=False, labelbottom=False)
    ax.tick_params(axis='y', labelsize=12, length=2.5, pad=2)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_linewidth(0.7)


def _draw_slot_occupancy(ax, matrix, method_color, n_cells, n_slots):
    """Plot the number of selected slots for each beam pattern."""
    slot_occupancy = np.count_nonzero(matrix, axis=1)
    ax.plot(slot_occupancy, np.arange(1, n_cells + 1),
            color='black', linewidth=1.25)
    ax.set_xlim(0, n_slots)
    ax.set_xticks([0, n_slots // 2, n_slots])
    ax.set_xlabel('Slot-Occupancy', fontsize=18,
                  fontname='Times New Roman', labelpad=3)
    ax.tick_params(axis='x', labelsize=12, length=2.5, pad=2)
    ax.tick_params(axis='y', left=False, labelleft=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.spines['left'].set_color('black')
    ax.spines['left'].set_linewidth(0.7)
    ax.spines['bottom'].set_color('black')
    ax.spines['bottom'].set_linewidth(0.7)


def _plot_beam_pattern_grid(method_specs, output_path, active_y_max):
    """Create one four-panel beam-pattern figure with shared visual scales."""
    fig = plt.figure(figsize=(18, 13))
    outer = fig.add_gridspec(2, 2, hspace=0.25, wspace=0.24)

    for idx, (method_name, matrix, method_color) in enumerate(method_specs):
        n_cells, n_slots = matrix.shape
        inner = outer[idx // 2, idx % 2].subgridspec(
            2, 2,
            height_ratios=[0.23, 1.0],
            width_ratios=[1.0, 0.17],
            hspace=0.055,
            wspace=0.045
        )
        ax_top = fig.add_subplot(inner[0, 0])
        ax_main = fig.add_subplot(inner[1, 0])
        ax_right = fig.add_subplot(inner[1, 1], sharey=ax_main)

        _draw_active_cell_profile(
            ax_top, matrix, method_color, n_slots, active_y_max
        )
        ax_top.text(
            0.0, 1.13, f'({chr(97 + idx)})  {method_name}',
            transform=ax_top.transAxes, ha='left', va='bottom',
            fontsize=18, fontname='Times New Roman'
        )
        _draw_beam_gantt(ax_main, matrix, method_color)
        _style_beam_pattern_axis(ax_main, n_slots, n_cells)
        _draw_slot_occupancy(
            ax_right, matrix, method_color, n_cells, n_slots
        )

    colorbar_ax = fig.add_axes([0.34, 0.960, 0.32, 0.016])
    connection_colorbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=_BEAM_CONNECTION_NORM,
                             cmap=_BEAM_CONNECTION_CMAP),
        cax=colorbar_ax, orientation='horizontal', ticks=np.arange(0, 9))
    connection_colorbar.ax.set_title('Beam Connection Count', fontsize=16, pad=7)
    connection_colorbar.ax.tick_params(labelsize=12, length=2.5, pad=2)
    connection_colorbar.outline.set_linewidth(0.7)
    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.07, top=0.91)
    fig.savefig(output_path, format='pdf', bbox_inches='tight')
    fig.savefig(output_path.replace('.pdf', '.svg'), format='svg', bbox_inches='tight')
    plt.close(fig)


all_beam_matrices = [
    connCount_C,
    connCount_C_ultimate,
    connCount_C_ultimate_MC,
    connCount_C_ultimate_TDM
]
active_cell_max = int(np.ceil(
    max(np.max(np.count_nonzero(matrix, axis=0)) for matrix in all_beam_matrices) / 5.0
) * 5)

output_path_figure3to6 = FIGURES_DIR + 'figure3to6.pdf'
_plot_beam_pattern_grid(
    [
        ('Coordinated BH', connCount_C, '#9467BD'),
        ('Proposed', connCount_C_ultimate, '#1F77B4'),
        ('Multi-Conn Heuristic', connCount_C_ultimate_MC, '#FF7F0E'),
        ('Terminal-Dominated Matching', connCount_C_ultimate_TDM, '#2CA02C')
    ],
    output_path_figure3to6,
    active_cell_max
)


mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['xtick.labelsize']=18
mpl.rcParams['ytick.labelsize']=18
def ecdf(data):
    x = np.sort(data)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y

# Figure 7: demand-dependent OCTR distributions.
df_rd = pd.read_csv(RESULTS_DIR + 'figure7.csv')

demand_levels = np.sort(df_rd['demand'].unique())
all_methods = df_rd['method'].unique()
methods_order = ['CPDA', 'Terminal-Dominated Matching', 'Multi-Conn Heuristic', 'Coordinated BH']
method_markers = {
    'CPDA': 'o',
    'Terminal-Dominated Matching': 'D',
    'Multi-Conn Heuristic': 's',
    'Coordinated BH': '^'
}



methods = [m for m in methods_order if m in all_methods]

COLORS3 = ['#1f77b4',  # blue
           '#2ca02c',  # green
           '#ff7f0e']  # orange

colors = COLORS3[:len(demand_levels)]




plt.figure(figsize=(8, 6))
for i, demand in enumerate(demand_levels):
    for method in methods:
        subset = df_rd[(df_rd['demand'] == demand) & (df_rd['method'] == method)]
        if subset.empty:
            continue
        x, y = ecdf(subset['value'].values)
        plt.plot(x, y, '-', color=colors[i], linewidth=1)

        sample_levels = np.arange(0.05, 1.05, 0.05)
        x_sample = np.interp(sample_levels, y, x)
        marker_style = method_markers.get(method, 'o')
        plt.plot(x_sample, sample_levels, marker_style,
                 markersize=6, markerfacecolor='none', color=colors[i], linewidth=1)

plt.xlabel('OCTR', fontsize=18)
plt.ylabel('CDF', fontsize=18)
plt.xlim([0, 1.2])
plt.grid(True)
plt.tight_layout()

from matplotlib.lines import Line2D

method_legend_elements = [
    Line2D([0], [0], marker=method_markers[m], linestyle='-',linewidth=1,
           color='k', markersize=6, markerfacecolor='none', label=m)
    for m in methods
]

demand_legend_elements = [
    Line2D([0], [0], color=colors[i], lw=1, label=f'{demand:.1f} Mbps')
    for i, demand in enumerate(demand_levels)
]
first_legend = plt.legend(handles=method_legend_elements, loc='lower right', fontsize=14)

plt.gca().add_artist(first_legend)
plt.legend(handles=demand_legend_elements, loc='upper left', fontsize=14)


output_path_figure7 = FIGURES_DIR + 'figure7.pdf'
plt.savefig(output_path_figure7, format='pdf', bbox_inches='tight')
plt.savefig(output_path_figure7.replace('.pdf', '.svg'), format='svg', bbox_inches='tight')



mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['xtick.labelsize']=18
mpl.rcParams['ytick.labelsize']=18

# Figure 8: demand-dependent Jain fairness distributions.
df_jain = pd.read_csv(RESULTS_DIR + 'figure8.csv')
df_jain['demand'] = pd.to_numeric(df_jain['demand'], errors='coerce')
demand_levels = np.sort(df_jain['demand'].dropna().unique())

plot_grouped_box_by_method(
    df=df_jain,
    group_col='demand',
    group_order=demand_levels,
    group_labels={d: f'{d:g} Mbps' for d in demand_levels},
    out_pdf=FIGURES_DIR + 'figure8.pdf',
    ylabel='Jain Fairness Index',
    xlabel='Traffic Demand Density',
    ylim=[0, 1.02],
    legend_loc='lower left'
)


mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18

def ecdf(data):
    data = np.asarray(data, dtype=float)
    data = data[~np.isnan(data)]
    x = np.sort(data)
    if len(x) == 0:
        return np.array([]), np.array([])
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y

methods_order = ['CPDA', 'Terminal-Dominated Matching', 'Multi-Conn Heuristic', 'Coordinated BH']
method_markers = {
    'CPDA': 'o',
    'Terminal-Dominated Matching': 'D',
    'Multi-Conn Heuristic': 's',
    'Coordinated BH': '^'
}

def color_for_config(cc):
    return 'g' if int(cc) == 1 else 'b'

def label_for_config(cc):
    return 'Phase-1 840' if int(cc) == 1 else 'Phase-2 2016'


# Figure 9: constellation-dependent OCTR distributions.
df_rd = pd.read_csv(RESULTS_DIR + 'figure9.csv')

df_rd['config'] = pd.to_numeric(df_rd['config'], errors='coerce')
df_rd = df_rd.dropna(subset=['config', 'value'])

configs = np.sort(df_rd['config'].unique())
present_methods = set(df_rd['method'].unique())
methods = [m for m in methods_order if m in present_methods]

plt.figure(figsize=(8, 6))

for cc in configs:
    for method in methods:
        subset = df_rd[(df_rd['config'] == cc) & (df_rd['method'] == method)]
        if subset.empty:
            continue

        x, y = ecdf(subset['value'].values)
        if len(x) == 0:
            continue

        plotColor = color_for_config(cc)

        plt.plot(x, y, '-', color=plotColor, linewidth=1)

        sample_levels = np.arange(0.05, 1.05, 0.05)
        x_sample = np.interp(sample_levels, y, x)
        mk = method_markers.get(method, 'o')
        plt.plot(
            x_sample, sample_levels, mk,
            markersize=6, markerfacecolor='none',
            color=plotColor, linewidth=1
        )

plt.xlabel('OCTR', fontsize=18)
plt.ylabel('CDF', fontsize=18)
plt.xlim([0, 1.1])
plt.grid(True)
plt.tight_layout()

method_legend_elements = [
    Line2D([0], [0], marker=method_markers[m], linestyle='-', linewidth=1,
           color='k', markersize=6, markerfacecolor='none', label=m)
    for m in methods
]

config_legend_elements = [
    Line2D([0], [0], color=color_for_config(cc), lw=1, label=label_for_config(cc))
    for cc in configs
]

first_legend = plt.legend(handles=method_legend_elements, loc='center right', fontsize=14)
plt.gca().add_artist(first_legend)
plt.legend(handles=config_legend_elements, loc='lower right', fontsize=14)

plt.savefig(FIGURES_DIR + 'figure9.pdf', format='pdf', bbox_inches='tight')
plt.savefig(FIGURES_DIR + 'figure9.svg', format='svg', bbox_inches='tight')
plt.close()


# Figure 10: constellation-dependent Jain fairness distributions.
df_jain = pd.read_csv(RESULTS_DIR + 'figure10.csv')
df_jain['config'] = pd.to_numeric(df_jain['config'], errors='coerce')
df_jain = df_jain.dropna(subset=['config', 'value'])

configs = np.sort(df_jain['config'].unique())
plot_grouped_box_by_method(
    df=df_jain,
    group_col='config',
    group_order=configs,
    group_labels={cc: label_for_config(cc) for cc in configs},
    out_pdf=FIGURES_DIR + 'figure10.pdf',
    ylabel='Jain Fairness Index',
    xlabel='Satellite Constellation',
    ylim=[0, 1.02],
    legend_loc='lower left'
)



mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18

def ecdf(data):
    data = np.asarray(data, dtype=float)
    data = data[~np.isnan(data)]
    x = np.sort(data)
    if len(x) == 0:
        return np.array([]), np.array([])
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y

# Figures 11-12: traffic-distribution comparisons.
scenario_order = ['Uniform', 'Moderate', 'Extreme']
scenario_colors = {
    'Uniform':  '#1f77b4',  # blue
    'Moderate': '#2ca02c',  # green
    'Extreme':  '#ff7f0e'   # orange
}

methods_order = ['CPDA', 'Terminal-Dominated Matching', 'Multi-Conn Heuristic', 'Coordinated BH']
method_markers = {
    'CPDA': 'o',
    'Terminal-Dominated Matching': 'D',
    'Multi-Conn Heuristic': 's',
    'Coordinated BH': '^'
}

def plot_cdf_from_longtable(csv_path, x_label, out_pdf,
                            legend1_loc='upper left', legend2_loc='lower right',
                            xlim=None):
    df = pd.read_csv(csv_path)

    df = df.dropna(subset=['scenario', 'method', 'value'])
    df['value'] = pd.to_numeric(df['value'], errors='coerce')
    df = df.dropna(subset=['value'])

    present_scenarios = set(df['scenario'].unique())
    present_methods = set(df['method'].unique())

    scenarios = [s for s in scenario_order if s in present_scenarios]
    methods = [m for m in methods_order if m in present_methods]

    plt.figure(figsize=(8, 6))

    for sc in scenarios:
        for method in methods:
            subset = df[(df['scenario'] == sc) & (df['method'] == method)]
            if subset.empty:
                continue

            x, y = ecdf(subset['value'].values)
            if len(x) == 0:
                continue

            c = scenario_colors.get(sc, '#1f77b4')
            mk = method_markers.get(method, 'o')

            plt.plot(x, y, '-', color=c, linewidth=1)

            sample_levels = np.arange(0.05, 1.05, 0.05)
            x_sample = np.interp(sample_levels, y, x)
            plt.plot(
                x_sample, sample_levels, mk,
                markersize=6, markerfacecolor='none',
                color=c, linewidth=1
            )

    plt.xlabel(x_label, fontsize=18)
    plt.ylabel('CDF', fontsize=18)
    if xlim is not None:
        plt.xlim(xlim)
    plt.grid(True)
    plt.tight_layout()

    method_legend_elements = [
        Line2D([0], [0],
               marker=method_markers[m],
               linestyle='-', linewidth=1,
               color='k', markersize=6,
               markerfacecolor='none',
               label=m)
        for m in methods
    ]

    scenario_legend_elements = [
        Line2D([0], [0], linestyle='-', color=scenario_colors[s], linewidth=1, label=s)
        for s in scenarios
    ]

    first_legend = plt.legend(handles=method_legend_elements, loc=legend1_loc, fontsize=14)
    plt.gca().add_artist(first_legend)
    plt.legend(handles=scenario_legend_elements, loc=legend2_loc, fontsize=14)

    plt.savefig(out_pdf, format='pdf', bbox_inches='tight')
    plt.savefig(out_pdf.replace('.pdf', '.svg'), format='svg', bbox_inches='tight')
    plt.close()


plot_cdf_from_longtable(
    csv_path=RESULTS_DIR + 'figure11.csv',
    x_label='OCTR',
    out_pdf=FIGURES_DIR + 'figure11.pdf',
    legend1_loc='center right',
    legend2_loc='lower right',
    xlim=[0, 1.1]
)

df_jain_scenario = pd.read_csv(RESULTS_DIR + 'figure12.csv')
plot_grouped_box_by_method(
    df=df_jain_scenario,
    group_col='scenario',
    group_order=scenario_order,
    group_labels={s: s for s in scenario_order},
    out_pdf=FIGURES_DIR + 'figure12.pdf',
    ylabel='Jain Fairness Index',
    xlabel='Traffic Demand Distribution',
    ylim=[0, 1.02],
    legend_loc='lower left'
)


import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18

# Per-cell rate, unmet demand, unused capacity, and system summary.
df_rate = pd.read_csv(RESULTS_DIR + 'figure_cell_rate.csv')
df_unmet = pd.read_csv(RESULTS_DIR + 'figure_cell_unmet.csv')
df_unused = pd.read_csv(RESULTS_DIR + 'figure_cell_unused.csv')
df_sum = pd.read_csv(RESULTS_DIR + 'figure_cell_summary.csv')

out_cell_all = FIGURES_DIR + 'figure_cell_rate_unmet_unused.pdf'

plot_cell_metrics_grouped_bars(df_rate, df_unmet, df_unused, out_cell_all)

print('[DONE] Saved:')
print(' -', out_cell_all)

print('\n[Summary Table] (SUM over setups, no averaging)')
print(df_sum)


import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18

# ASE comparison with and without fragmentation.
csv_path = RESULTS_DIR + 'figure_ASE_three_methods.csv'
df = pd.read_csv(csv_path)

df['method'] = df['method'].astype(str).str.strip()

df['method'] = df['method'].replace({
    'CPDA (NoFrag)': 'CPDA (Without Fragmentation)',
    'CPDA(NoFrag)': 'CPDA (Without Fragmentation)',
    'CPDA (No Fragmentation)': 'CPDA (Without Fragmentation)',
    'CPDA(No Fragmentation)': 'CPDA (Without Fragmentation)',
    'CPDA (Without Fragmentation)': 'CPDA (Without Fragmentation)',
})

methods_order = ['CPDA', 'CPDA (Without Fragmentation)', 'Coordinated BH']
ase_colors = {
    'CPDA': METHOD_COLORS['CPDA'],
    'CPDA (Without Fragmentation)': METHOD_COLORS['CPDA (Without Fragmentation)'],
    'Coordinated BH': METHOD_COLORS['Coordinated BH']
}
ase_linestyles = {
    'CPDA': '-',
    'CPDA (Without Fragmentation)': '--',
    'Coordinated BH': '-'
}

out_pdf = FIGURES_DIR + 'figure_ASE_three_methods.pdf'
plot_box_by_method(
    df=df,
    value_col='ASE',
    out_pdf=out_pdf,
    ylabel='Area Spectral Efficiency (bit/s/Hz/km$^2$)',
    methods_order=methods_order,
    colors=ase_colors,
    linestyles=ase_linestyles
)
print('[DONE] Saved:', out_pdf)

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18

# Bandwidth-dependent ASE gap.
csv_path = RESULTS_DIR + 'figure_ASE_bandwidth_gap.csv'
df = pd.read_csv(csv_path)

base_out = FIGURES_DIR
out_pdf = base_out + 'figure_ASE_bandwidth_gap.pdf'
plot_bandwidth_ase_boxes(df, out_pdf)
print('[DONE] Saved:', out_pdf)



# Figures 17-20: layout and scheduling runtime distributions.
runtime17_df = pd.read_csv(RESULTS_DIR + 'figure17_runtime.csv')
runtime18_df = pd.read_csv(RESULTS_DIR + 'figure18_runtime.csv')

runtime17_pdf = FIGURES_DIR + 'figure17_runtime.pdf'
plot_runtime_boxplot(runtime17_df, runtime17_pdf)
print('[DONE] Saved:', runtime17_pdf)

runtime18_pdf = FIGURES_DIR + 'figure18_runtime.pdf'
plot_runtime_boxplot(runtime18_df, runtime18_pdf)
print('[DONE] Saved:', runtime18_pdf)

comparison_order = ['CPDA', 'Layout Generation']
runtime19_pdf = FIGURES_DIR + 'figure19_runtime.pdf'
plot_runtime_boxplot(runtime17_df, runtime19_pdf, method_order=comparison_order)
print('[DONE] Saved:', runtime19_pdf)

runtime20_pdf = FIGURES_DIR + 'figure20_runtime.pdf'
plot_runtime_boxplot(runtime18_df, runtime20_pdf, method_order=comparison_order)
print('[DONE] Saved:', runtime20_pdf)
