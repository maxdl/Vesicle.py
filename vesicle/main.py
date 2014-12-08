# -*- coding: utf-8 -*-

import sys
import os.path
import time
import core
import geometry
import file_io
import version
import stringconv as sc

#
# Functions
#


def save_output(profileli, opt):
    """ Save a summary of results of evaluated profiles
    """
    def m(x, pixelwidth):
        return geometry.to_metric_units(x, pixelwidth)

    def m2(x, pixelwidth):
        # For area units
        return geometry.to_metric_units(x, pixelwidth**2)

    def na(x):
        if x in (None, -1):
            return "N/A"
        else:
            return x

    def write_session_summary():
        with file_io.FileWriter("session.summary", opt) as f:
            f.writerow(["%s version:" % version.title,
                       "%s (Last modified %s %s, %s)"
                       % ((version.version,) + version.date)])
            f.writerow(["Number of evaluated profiles:", len(eval_proli)])
            if err_fli:
                f.writerow(["Number of non-evaluated profiles:",
                            len(err_fli)])
            f.writerow(["Metric unit:", eval_proli[0].metric_unit])
            f.writerow(["Allow vesicle overlap:",
                        sc.yes_or_no(opt.allow_vesicle_overlap)])
            f.writerow(["Prioritize vesicle lumen over membrane:",
                        sc.yes_or_no(opt.prioritize_lumen)])
            f.writerow(["Intervesicle distances calculated:",
                        sc.yes_or_no(opt.determine_intervesicle_dists)])
            if opt.determine_intervesicle_dists:
                f.writerow(["Intervesicle distance mode:",
                            opt.intervesicle_dist_mode])
                f.writerow(["Shortest intervesicle distances:",
                            sc.yes_or_no(opt.intervesicle_shortest_dist)])
                f.writerow(["Lateral intervesicle distances:",
                            sc.yes_or_no(opt.intervesicle_lateral_dist)])
            f.writerow(["Spatial resolution:", opt.spatial_resolution,
                        eval_proli[0].metric_unit])
            f.writerow(["Shell width:", opt.shell_width,
                        eval_proli[0].metric_unit])
            f.writerow(["Interpoint distances calculated:",
                        sc.yes_or_no(opt.determine_interpoint_dists)])
            if opt.determine_interpoint_dists:
                f.writerow(["Interpoint distance mode:",
                            opt.interpoint_dist_mode])
                f.writerow(["Shortest interpoint distances:",
                            sc.yes_or_no(opt.interpoint_shortest_dist)])
                f.writerow(["Lateral interpoint distances:",
                            sc.yes_or_no(opt.interpoint_lateral_dist)])
            if clean_fli:
                f.writerow(["Input files processed cleanly:"])
                f.writerows([[fn] for fn in clean_fli])
            if nop_fli:
                f.writerow(["Input files processed but which generated no "
                            "point distances:"])
                f.writerows([[fn] for fn in nop_fli])
            if warn_fli:
                f.writerow(["Input files processed but which generated "
                            "warnings (see log for details):"])
                f.writerows([[fn] for fn in warn_fli])
            if err_fli:
                f.writerow(["Input files not processed or not included in "
                            "summary (see log for details):"])
                f.writerows([[fn] for fn in err_fli])

    def write_profile_summary():
        table = [["Perimeter",
                  "Area",
                  "Vesicles",
                  "Vesicle density",
                  "Vesicle perimeter (total)",
                  "Vesicle area (total)"]]
        if opt.points_found:
            table[0].extend(["Points (total)",
                             "Points within profile",
                             "Area density of points within profile",
                             "Vesicle-associated points"])
        if opt.use_random:
            table[0].extend(["Random points (total)",
                             "Random points in profile",
                             "Vesicle-associated random points"])
        table[0].extend(["Profile ID",
                         "Input file",
                         "Comment"])
        for pro in eval_proli:
            row = [m(pro.path.perimeter(), pro.pixelwidth),
                   m2(pro.path.area(), pro.pixelwidth),
                   len(pro.vli),
                   1e6 * (len([v for v in pro.vli])
                          / m2(pro.path.area(), pro.pixelwidth)),
                   m(sum([v.perimeter() for v in pro.vli]),
                     pro.pixelwidth),
                   m2(sum([v.area() for v in pro.vli]), pro.pixelwidth)]
            if opt.points_found:
                row.extend([len(pro.pli),
                           len([p for p in pro.pli if p.is_within_profile]),
                           1e6*(len([p for p in pro.pli if p.is_within_profile])
                                / m2(pro.path.area(), pro.pixelwidth)),
                           len([p for p in pro.pli
                                if p.is_associated_with_vesicle])])
            if opt.use_random:
                row.extend([len(pro.randomli),
                            len([r for r in pro.randomli
                                if r.is_within_profile]),
                            len([r for r in pro.randomli
                                if r.is_associated_with_vesicle])])
            row.extend([pro.ID,
                        os.path.basename(pro.inputfn),
                        pro.comment])
            table.append(row)
        with file_io.FileWriter("profile.summary", opt) as f:
            f.writerows(table)

    def write_vesicle_summary():
        with file_io.FileWriter("vesicle.summary", opt) as f:
            f.writerow(["Vesicle number (as appearing in input file)",
                        "Perimeter",
                        "Area",
                        "Diameter of circle of equal perimeter",
                        "Diameter of circle of equal area",
                        "Circularity",
                        "Convexity",
                        "Profile ID",
                        "Input file",
                        "Comment"])
            f.writerows([[n+1,
                          m(v.perimeter(), pro.pixelwidth),
                          m2(v.area(), pro.pixelwidth),
                          m(v.diameter_of_equal_perimeter_circle(),
                            pro.pixelwidth),
                          m(v.diameter_of_equal_area_circle(),
                            pro.pixelwidth),
                          v.circularity(),
                          v.convexity(),
                          pro.ID,
                          os.path.basename(pro.inputfn),
                          pro.comment]
                        for pro in eval_proli
                        for n, v in enumerate(pro.vli)])

    def write_point_summary(ptype):
        if ptype == "point":
            if not opt.points_found:
                return
            pli = "pli"
            pstr = "point"
        elif ptype == "random":
            if not opt.use_random:
                return
            pli = "randomli"
            pstr = "point"
        else:
            return
        with file_io.FileWriter("%s.summary" % ptype, opt) as f:
            f.writerow(["%s number (as appearing in input file)"
                        % pstr.capitalize(),
                        "Distance to nearest vesicle border",
                        "Distance to nearest vesicle center",
                        "Within vesicle",
                        "Vesicle-associated",
                        "Distance to profile border",                                                
                        "Within profile",
                        "Profile border-associated",
                        "Profile-associated",                      
                        "Profile ID",
                        "Input file",
                        "Comment"])
            f.writerows([[n+1,
                          m(p.dist_to_nearest_vesicle_border, pro.pixelwidth),
                          m(p.dist_to_nearest_vesicle_center, pro.pixelwidth),
                          sc.yes_or_no(p.is_within_vesicle),
                          sc.yes_or_no(p.is_associated_with_vesicle),
                          m(p.dist_to_path, pro.pixelwidth),                                                    
                          sc.yes_or_no(p.is_within_profile),
                          sc.yes_or_no(p.is_associated_with_path),
                          sc.yes_or_no(p.is_within_profile or
                                       p.is_associated_with_path),
                          pro.ID,
                          os.path.basename(pro.inputfn),
                          pro.comment] for pro in eval_proli for n, p in
                         enumerate(pro.__dict__[pli])])

    def write_inter_summaries(comp):

        def _m(x):
            return m(x, pro.pixelwidth)

        if not opt.__dict__["determine_inter%s_dists" % comp]:
            return
        original_inter_relations = opt.__dict__["inter%s_relations" % comp]
        dist_mode = opt.__dict__["inter%s_dist_mode" % comp]
        shortest_dist = opt.__dict__["inter%s_shortest_dist" % comp]
        lateral_dist = opt.__dict__["inter%s_lateral_dist" % comp]
        inter_rels = dict([(key, val)
                          for key, val in original_inter_relations.items()
                          if val])
        if not opt.use_random:
            for key in original_inter_relations.keys():
                if "random" in key:
                    del inter_rels[key]
        if len(inter_rels) == 0 or not (shortest_dist or lateral_dist):
            return
        table = []
        if dist_mode == 'all':
            s = "all distances"
        else:
            s = "nearest neighbour distances"
        table.append(["Mode: " + s])
        keyli = inter_rels.keys()
        headerli = [k for k in keyli]
        if shortest_dist and lateral_dist:
            headerli.extend(headerli)
        topheaderli = []
        if shortest_dist:
            topheaderli.append("Shortest distances")
            if lateral_dist:
                topheaderli.extend([""] * (len(inter_rels) - 1))
        if lateral_dist:
            topheaderli.append("Lateral distances along profile border "
                               "membrane")
        table.extend([topheaderli, headerli])
        cols = [[] for c in headerli]
        lateral_col0 = len(keyli) if shortest_dist else 0
        for pro in eval_proli:
            for n, di in enumerate([pro.interdistlis[key] for key in keyli]):
                if shortest_dist:
                    cols[n].extend(map(_m, di["shortest"]))
                if lateral_dist:
                    cols[n + lateral_col0].extend(map(_m, di["lateral"]))
        # transpose cols and append to table
        table.extend(map(lambda *col: [e if e is not None else "" for e in col],
                         *cols))
        with file_io.FileWriter("inter%s.summary" % comp, opt) as f:
            f.writerows(table)

    sys.stdout.write("\nSaving summaries...\n")
    opt.save_result = {'any_saved': False, 'any_err': False}
    eval_proli = [profile for profile in profileli if not profile.errflag]
    clean_fli = [profile.inputfn for profile in profileli
                 if not (profile.errflag or profile.warnflag)]
    warn_fli = [profile.inputfn for profile in profileli if profile.warnflag]
    err_fli = [profile.inputfn for profile in profileli if profile.errflag]
    nop_fli = [profile.inputfn for profile in profileli if not profile.pli]
    for profile in profileli:
        if len(profile.pli) > 0:
            opt.points_found = True
            break
    else:
        opt.points_found = False
    write_session_summary()
    write_profile_summary()
    write_vesicle_summary()
    write_point_summary("point")
    write_point_summary("random")
    write_inter_summaries("vesicle")
    write_inter_summaries("point")
    if opt.save_result['any_err']:
        sys.stdout.write("Note: One or more summaries could not be saved.\n")
    if opt.save_result['any_saved']:
        sys.stdout.write("Done.\n")
    else:
        sys.stdout.write("No summaries saved.\n")


def reset_options(opt):
    """ Deletes certain options that should always be set anew for each run
        (each time the "Start" button is pressed)
    """
    opt.metric_unit = None
    opt.use_random = None


def show_options(opt):
    sys.stdout.write("{} version: {} (Last modified {} {}, {})\n".format(
                     version.title, version.version, *version.date))
    sys.stdout.write("Output file format: %s\n" % opt.output_file_format)
    sys.stdout.write("Suffix of output files: %s\n"
                     % opt.output_filename_suffix)
    sys.stdout.write("Output directory: %s\n" % opt.output_dir)
    sys.stdout.write("Allow vesicle overlap: %s\n"
                     % sc.yes_or_no(opt.allow_vesicle_overlap))
    sys.stdout.write("Prioritize vesicle lumen over membrane: %s\n"
                     % sc.yes_or_no(opt.prioritize_lumen))
    sys.stdout.write("Spatial resolution: %d\n" % opt.spatial_resolution)
    sys.stdout.write("Shell width: %d metric units\n" % opt.shell_width)
    sys.stdout.write("Intervesicle distances calculated: %s\n"
                     % sc.yes_or_no(opt.determine_intervesicle_dists))
    if opt.determine_intervesicle_dists:
        sys.stdout.write("Intervesicle distance mode: %s\n"
                         % opt.intervesicle_dist_mode.capitalize())
        sys.stdout.write("Shortest intervesicle distances: %s\n"
                         % sc.yes_or_no(opt.intervesicle_shortest_dist))
        sys.stdout.write("Lateral intervesicle distances: %s\n"
                         % sc.yes_or_no(opt.intervesicle_lateral_dist))    
    sys.stdout.write("Interpoint distances calculated: %s\n"
                     % sc.yes_or_no(opt.determine_interpoint_dists))
    if opt.determine_interpoint_dists:
        sys.stdout.write("Interpoint distance mode: %s\n"
                         % opt.interpoint_dist_mode.capitalize())
        sys.stdout.write("Shortest interpoint distances: %s\n"
                         % sc.yes_or_no(opt.interpoint_shortest_dist))
        sys.stdout.write("Lateral interpoint distances: %s\n"
                         % sc.yes_or_no(opt.interpoint_lateral_dist))


def get_output_format(opt):
    if opt.output_file_format == 'excel':
        import imp
        try:
            imp.find_module("pyExcelerator")
        except ImportError:
            sys.stdout.write("Unable to write Excel files: resorting to csv "
                             "format.\n")
            opt.output_file_format = "csv"
    if opt.output_file_format == 'csv':
        opt.output_filename_ext = ".csv"
        opt.csv_format = {'dialect': 'excel', 'lineterminator': '\n',
                          'encoding': sys.getfilesystemencoding()}
        if opt.csv_delimiter == 'tab':
            opt.csv_format['delimiter'] = '\t'
    if opt.output_filename_date_suffix:
        from datetime import date
        opt.output_filename_suffix = "." + date.today().isoformat()
    if opt.output_filename_other_suffix != '':
        opt.output_filename_suffix += "." + opt.output_filename_other_suffix

      
def main_proc(parent):
    """ Process profile data files
    """
    opt = parent.opt
    if not opt.input_file_list:
        sys.stdout.write("No input files.\n")
        return 0
    i, n = 0, 0
    profileli = []
    sys.stdout.write("--- Session started %s local time ---\n" % time.ctime())
    # Remove duplicate filenames
    for f in opt.input_file_list:
        if opt.input_file_list.count(f) > 1:
            sys.stdout.write("Duplicate input filename %s:\n   => "
                             "removing first occurrence in list\n" % f)
            opt.input_file_list.remove(f)
    get_output_format(opt)
    reset_options(opt)
    show_options(opt)
    while True:
        if i < len(opt.input_file_list):
            inputfn = opt.input_file_list[i]
            i += 1
        else: 
            sys.stdout.write("\nNo more input files...\n")
            break
        parent.process_queue.put(("new_file", inputfn))
        profileli.append(core.ProfileData(inputfn, opt))
        profileli[-1].process()
        if opt.stop_requested:
            sys.stdout.write("\n--- Session aborted by user %s local time ---\n"
                             % time.ctime())
            return 3
        if not profileli[-1].errflag:
            n += 1
            if profileli[-1].warnflag:
                sys.stdout.write("Warning(s) found while processing "
                                 "input file.\n")
                continue
        else:
            sys.stdout.write("Error(s) found while processing input file =>\n"
                             "  => No distances could be determined.\n")
            continue
    # no more input files
    errfli = [pro.inputfn for pro in profileli if pro.errflag]
    warnfli = [pro.inputfn for pro in profileli if pro.warnflag]
    if errfli:
        sys.stdout.write("\n%s input %s generated one or more errors:\n"
                         % (sc.plurality("This", len(errfli)),
                            sc.plurality("file", len(errfli))))
        sys.stdout.write("%s\n" % "\n".join([fn for fn in errfli]))
    if warnfli:
        sys.stdout.write("\n%s input %s generated one or more warnings:\n"
                         % (sc.plurality("This", len(warnfli)),
                            sc.plurality("file", len(warnfli))))
        sys.stdout.write("%s\n" % "\n".join([fn for fn in warnfli]))
    if n > 0:
        parent.process_queue.put(("saving_summaries", ""))
        save_output(profileli, opt)
    else:
        sys.stdout.write("\nNo files processed.\n")
    sys.stdout.write("--- Session ended %s local time ---\n" % time.ctime())
    parent.process_queue.put(("done", ""))
    opt.reset()
    if errfli: 
        return 0
    elif warnfli: 
        return 2
    else: 
        return 1
# End of main.py
