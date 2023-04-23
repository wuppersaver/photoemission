import math
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np


#Everything below this is to plot bandstructure plots with functions adapted from the pymatgen module

def get_bs_plot(
        bandstructure,
        zero_to_efermi=True,
        ylim=None,
        smooth=False,
        vbm_cbm_marker=False,
        smooth_tol=0,
        smooth_k=3,
        smooth_np=100,
        bs_labels=[],
        size = [10,8]
    ):
        plt.style.use('seaborn-poster')
        plt.rcParams["font.family"] = "serif"
        pretty = pretty_plot(size[0],size[1])
        bs_array = [bandstructure]

        if isinstance(smooth, bool):
            smooth = [smooth] * len(bs_array)

        handles = []
        vbm_min, cbm_max = [], []

        colors = list(pretty.rcParams["axes.prop_cycle"].by_key().values())[0]
        for ibs, bs in enumerate(bs_array):

            # set first bs in the list as ref for rescaling the distances of the other bands
            bs_ref = bs_array[0] if len(bs_array) > 1 and ibs > 0 else None

            if smooth[ibs]:
                # interpolation works good on short segments like branches
                data = bs_plot_data(zero_to_efermi, bs, bs_ref, split_branches=True)
            else:
                data = bs_plot_data(zero_to_efermi, bs, bs_ref, split_branches=False)

            # remember if one bs is a metal for setting the ylim later
            one_is_metal = False
            if not one_is_metal and data["is_metal"]:
                one_is_metal = data["is_metal"]

            # remember all the cbm and vbm for setting the ylim later
            if not data["is_metal"]:
                cbm_max.append(data["cbm"][0][1])
                vbm_min.append(data["vbm"][0][1])

            for sp in bs.bands.keys():
                ls = "-" if str(sp) == "1" else "--"

                if bs_labels != []:
                    bs_label = f"{bs_labels[ibs]} {sp.name}"
                else:
                    bs_label = f"Band {ibs} {sp.name}"

                handles.append(mlines.Line2D([], [], lw=0.2, ls=ls, color=colors[ibs], label=bs_label))

                distances, energies = data["distances"], data["energy"][str(sp)]

                for dist, ene in zip(distances, energies):
                    pretty.plot(dist, ene.T, ls=ls)

            # plot markers for vbm and cbm
            if vbm_cbm_marker:
                for cbm in data["cbm"]:
                    pretty.scatter(cbm[0], cbm[1], color="r", marker="o", s=100)
                for vbm in data["vbm"]:
                    pretty.scatter(vbm[0], vbm[1], color="g", marker="o", s=100)

            # Draw Fermi energy, only if not the zero
            if not zero_to_efermi:
                ef = bs.efermi
                pretty.axhline(ef, lw=1, ls="-.", color=colors[ibs])

        # defaults for ylim
        e_min = -4
        e_max = 4
        if one_is_metal:
            e_min = -10
            e_max = 10

        if ylim is None:
            if zero_to_efermi:
                if one_is_metal:
                    # Plot A Metal
                    pretty.ylim(e_min, e_max)
                else:
                    pretty.ylim(e_min, max(cbm_max) + e_max)
            else:
                all_efermi = [b.efermi for b in bs_array]
                ll = min([min(vbm_min), min(all_efermi)])
                hh = max([max(cbm_max), max(all_efermi)])
                pretty.ylim(ll + e_min, hh + e_max)
        else:
            pretty.ylim(ylim)
        pretty = maketicks(bs, pretty)

        # Main X and Y Labels
        #plt.xlabel(r"$\mathrm{Wave\ Vector}$", fontsize=30)
        pretty.xlabel('')
        ylabel = r"$\mathrm{E\ -\ E_f\ (eV)}$" if zero_to_efermi else r"$\mathrm{Energy\ (eV)}$"
        pretty.ylabel(ylabel, fontsize=15)

        # X range (K)
        # last distance point
        x_max = data["distances"][-1][-1]
        pretty.xlim(0, x_max)

        #plt.legend(handles=handles)

        pretty.tight_layout()

        # auto tight_layout when resizing or pressing t
        def fix_layout(event):
            if (event.name == "key_press_event" and event.key == "t") or event.name == "resize_event":
                pretty.gcf().tight_layout()
                pretty.gcf().canvas.draw()

        pretty.gcf().canvas.mpl_connect("key_press_event", fix_layout)
        pretty.gcf().canvas.mpl_connect("resize_event", fix_layout)

        return pretty

def pretty_plot(width=8, height=None, plt=None, dpi=300, color_cycle=("qualitative", "Set1_9")):
    """
    Provides a publication quality plot, with nice defaults for font sizes etc.

    Args:
        width (float): Width of plot in inches. Defaults to 8in.
        height (float): Height of plot in inches. Defaults to width * golden
            ratio.
        plt (matplotlib.pyplot): If plt is supplied, changes will be made to an
            existing plot. Otherwise, a new plot will be created.
        dpi (int): Sets dot per inch for figure. Defaults to 300.
        color_cycle (tuple): Set the color cycle for new plots to one of the
            color sets in palettable. Defaults to a qualitative Set1_9.

    Returns:
        Matplotlib plot object with properly sized fonts.
    """
    ticksize = int(width)*1.6

    golden_ratio = (math.sqrt(5) - 1) / 2

    if not height:
        height = int(width * golden_ratio)

    if plt is None:
        import importlib

        import matplotlib.pyplot as plt

        mod = importlib.import_module(f"palettable.colorbrewer.{color_cycle[0]}")
        colors = getattr(mod, color_cycle[1]).mpl_colors
        from cycler import cycler

        plt.figure(figsize=(width, height), facecolor="w", dpi=dpi)
        ax = plt.gca()
        ax.set_prop_cycle(cycler("color", colors))
    else:
        fig = plt.gcf()
        fig.set_size_inches(width, height)
    plt.xticks(fontsize=ticksize)
    plt.yticks(fontsize=ticksize)

    ax = plt.gca()
    ax.set_title(ax.get_title(), size=width * 1.6)

    labelsize = int(width*4)

    ax.set_xlabel(ax.get_xlabel(), size=labelsize)
    ax.set_ylabel(ax.get_ylabel(), size=labelsize)

    return plt

def bs_plot_data(zero_to_efermi=True,bs = None, bs_ref=None, split_branches=True):
        """
        Get the data nicely formatted for a plot

        Args:
            zero_to_efermi: Automatically subtract off the Fermi energy from the
                eigenvalues and plot.
            bs: the bandstructure to get the data from. If not provided, the first
                one in the self._bs list will be used.
            bs_ref: is the bandstructure of reference when a rescale of the distances
                is need to plot multiple bands
            split_branches: if True distances and energies are split according to the
                branches. If False distances and energies are split only where branches
                are discontinuous (reducing the number of lines to plot).

        Returns:
            dict: A dictionary of the following format:
            ticks: A dict with the 'distances' at which there is a kpoint (the
            x axis) and the labels (None if no label).
            energy: A dict storing bands for spin up and spin down data
            {Spin:[np.array(nb_bands,kpoints),...]} as a list of discontinuous kpath
            of energies. The energy of multiple continuous branches are stored together.
            vbm: A list of tuples (distance,energy) marking the vbms. The
            energies are shifted with respect to the fermi level is the
            option has been selected.
            cbm: A list of tuples (distance,energy) marking the cbms. The
            energies are shifted with respect to the fermi level is the
            option has been selected.
            lattice: The reciprocal lattice.
            zero_energy: This is the energy used as zero for the plot.
            band_gap:A string indicating the band gap and its nature (empty if
            it's a metal).
            is_metal: True if the band structure is metallic (i.e., there is at
            least one band crossing the fermi level).
        """

        energies = {str(sp): [] for sp in bs.bands.keys()}

        bs_is_metal = bs.is_metal()

        if not bs_is_metal:
            vbm = bs.get_vbm()
            cbm = bs.get_cbm()

        zero_energy = 0.0
        if zero_to_efermi:
            if bs_is_metal:
                zero_energy = bs.efermi
            else:
                zero_energy = vbm["energy"]

        # rescale distances when a bs_ref is given as reference,
        # and when bs and bs_ref have different points in branches.
        # Usually bs_ref is the first one in self._bs list is bs_ref
        distances = bs.distance
        if bs_ref is not None:
            if bs_ref.branches != bs.branches:
                distances = bs._rescale_distances(bs_ref, bs)

        if split_branches:
            steps = [br["end_index"] + 1 for br in bs.branches][:-1]
        else:
            # join all the continuous branches
            # to reduce the total number of branches to plot
            steps = get_branch_steps(bs.branches)[1:-1]

        distances = np.split(distances, steps)
        for sp in bs.bands.keys():
            energies[str(sp)] = np.hsplit(bs.bands[sp] - zero_energy, steps)

        ticks = get_ticks(bs)

        vbm_plot = []
        cbm_plot = []
        bg_str = ""

        if not bs_is_metal:
            for index in cbm["kpoint_index"]:
                cbm_plot.append(
                    (
                        bs.distance[index],
                        cbm["energy"] - zero_energy if zero_to_efermi else cbm["energy"],
                    )
                )

            for index in vbm["kpoint_index"]:
                vbm_plot.append(
                    (
                        bs.distance[index],
                        vbm["energy"] - zero_energy if zero_to_efermi else vbm["energy"],
                    )
                )

            bg = bs.get_band_gap()
            direct = "Indirect"
            if bg["direct"]:
                direct = "Direct"

            bg_str = f"{direct} {bg['transition']} bandgap = {bg['energy']}"

        return {
            "ticks": ticks,
            "distances": distances,
            "energy": energies,
            "vbm": vbm_plot,
            "cbm": cbm_plot,
            "lattice": bs.lattice_rec.as_dict(),
            "zero_energy": zero_energy,
            "is_metal": bs_is_metal,
            "band_gap": bg_str,
        }

def get_branch_steps(branches):
        """
        Method to find discontinuous branches
        """
        steps = [0]
        for b1, b2 in zip(branches[:-1], branches[1:]):
            if b2["name"].split("-")[0] != b1["name"].split("-")[-1]:
                steps.append(b2["start_index"])
        steps.append(branches[-1]["end_index"] + 1)
        return steps

def get_ticks(bs):
        """
        Get all ticks and labels for a band structure plot.

        Returns:
            dict: A dictionary with 'distance': a list of distance at which
            ticks should be set and 'label': a list of label for each of those
            ticks.
        """
        bs = [bs]
        bs = bs[0] if isinstance(bs, list) else bs
        ticks, distance = [], []
        #print(bs.branches)
        for br in bs.branches:
            s, e = br["start_index"], br["end_index"]

            labels = br["name"].split("-")

            # skip those branches with only one point
            if labels[0] == labels[1]:
                continue

            # add latex $$
            for i, l in enumerate(labels):
                if l.startswith("\\") or "_" in l:
                    labels[i] = "$" + l + "$"

            # If next branch is not continuous,
            # join the first lbl to the previous tick label
            # and add the second lbl to ticks list
            # otherwise add to ticks list both new labels.
            # Similar for distances.
            if ticks and labels[0] != ticks[-1]:
                ticks[-1] += "$\\mid$" + labels[0]
                ticks.append(labels[1])
                distance.append(bs.distance[e])
            else:
                ticks.extend(labels)
                distance.extend([bs.distance[s], bs.distance[e]])

        return {"distance": distance, "label": ticks}

def maketicks(bs, plt):
    """
    utility private method to add ticks to a band structure
    """
    ticks = get_ticks(bs)
    print(ticks)
    # Sanitize only plot the uniq values
    uniq_d = []
    uniq_l = []
    temp_ticks = list(zip(ticks["distance"], ticks["label"]))
    for i, t in enumerate(temp_ticks):
        if i == 0:
            uniq_d.append(t[0])
            uniq_l.append(t[1])
            
        else:
            if t != temp_ticks[i-1]:
                uniq_d.append(t[0])
                uniq_l.append(t[1])

    plt.gca().set_xticks(uniq_d)
    plt.gca().set_xticklabels(uniq_l)

    for i in range(len(ticks["label"])):
        if ticks["label"][i] is not None:
            # don't print the same label twice
            if i != 0:
                plt.axvline(ticks["distance"][i], color="k",lw = 1)
            else:
                print(ticks["label"][i])
                if ticks["label"][i] != ticks["label"][i]:
                    plt.axvline(ticks["distance"][i], color="k",lw = 1)
    return plt