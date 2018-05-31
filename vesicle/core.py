import math
import sys
from . import geometry
from . import file_io


# Convenience functions

def dot_progress(line_length=80, char='.', reset=False):
    """Simple progress indicator on sys.stdout"""
    if not hasattr(dot_progress, 'counter'):
        dot_progress.counter = 0
    if reset:
        dot_progress.counter = 0
        sys.stdout.write('\n')
    dot_progress.counter += 1
    sys.stdout.write(char)
    if dot_progress.counter == line_length:
        dot_progress.counter = 0
        sys.stdout.write('\n')


def lazy_property(fn):
    """Decorator that makes a property lazily evaluated.
       From https://stevenloria.com/lazy-properties/.
    """
    attr_name = '_lazy_' + fn.__name__

    @property
    def _lazy_property(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazy_property


#
# Classes
#


class Point(geometry.Point):
    def __init__(self, x=None, y=None, ptype='', profile=None):
        if isinstance(x, geometry.Point):
            geometry.Point.__init__(self, x.x, x.y)
        else:
            geometry.Point.__init__(self, x, y)
        self.profile = profile
        if self.profile is not None:
            self.opt = self.profile.opt
        else:
            self.opt = None
        self.discard = False
        self.ptype = ptype
        self.nearest_neighbour_dist = None
        self.nearest_neighbour_point = geometry.Point()
        self.nearest_lateral_neighbour_dist = None
        self.nearest_lateral_neighbour_point = geometry.Point()
        self.nearest_neighbour = geometry.Point()

    def determine_stuff(self):
        """Determine general stuff for a point, including distance to path.
         Also mark the point for discarding if it is not valid.
        """

        def mark_to_discard(msg):
            if self.ptype == 'particle':  # don't warn if random points
                profile_message("Discarding particle at %s: %s" % (self, msg))
            self.discard = True
            self.profile.n_discarded[self.ptype] += 1
            return

        if self.is_within_hole:
            mark_to_discard("Located within a profile hole")
            return
        if not (self.is_within_profile or self.is_within_shell):
            mark_to_discard("Located outside the shell")
            return
        # This is to force the computation of this lazy property here
        __ = self.is_associated_with_path
        __ = self.dist_to_nearest_vesicle_center
        __ = self.dist_to_nearest_vesicle_border
        __ = self.is_within_vesicle

    @lazy_property
    def dist_to_path(self):
        """Return distance to profile border"""
        _dist_to_path = self.perpend_dist_closed_path(self.profile.path)
        if not self.is_within_profile:
            _dist_to_path = -_dist_to_path
        return _dist_to_path

    @lazy_property
    def is_within_hole(self):
        """Determine if self is inside a profile hole."""
        within_hole = False
        for h in self.profile.holeli:
            if self.is_within_polygon(h):
                within_hole = True
            else:
                within_hole = False
        return within_hole

    @lazy_property
    def is_within_profile(self):
        """Determine if self is inside profile, excluding holes."""
        if self.is_within_polygon(self.profile.path) and not self.is_within_hole:
            return True
        else:
            return False

    @lazy_property
    def is_within_shell(self):
        """Determine whether self is within shell"""
        return (not self.is_within_profile and
                abs(self.dist_to_path) < geometry.to_pixel_units(
                    self.opt.shell_width, self.profile.pixelwidth))

    @lazy_property
    def is_associated_with_path(self):
        """Determine whether self is associated with the profile
        border, i e, is within a distance of it that is less than
        the spatial resolution"""
        if (abs(self.dist_to_path) <= geometry.to_pixel_units(
                self.opt.spatial_resolution, self.profile.pixelwidth)):
            return True
        else:
            return False

    @lazy_property
    def is_associated_with_profile(self):
        """Determine whether self is within the profile or
        associated with the profile border"""
        if self.is_within_profile or self.is_associated_with_path:
            return True
        else:
            return False

    @lazy_property
    def is_within_vesicle(self):
        """Return True if self is within a vesicle"""
        for v in self.profile.vli:
            if self.is_within_polygon(v):
                return True
        return False

    @lazy_property
    def is_associated_with_vesicle(self):
        """Return True if self is within or associated with the border of a vesicle"""
        if  self.dist_to_nearest_vesicle_border  <= geometry.to_pixel_units(
            self.profile.opt.spatial_resolution, self.profile.pixelwidth):
            return True
        else:
            return False

    @lazy_property
    def dist_to_nearest_vesicle_center(self):
        """Return distance to nearest vesicle center"""
        vli = []
        if self.opt.prioritize_lumen:
            vli = [v for v in self.profile.vli if self.is_within_polygon(v)]
        if not vli:
            vli = self.profile.vli
        mindist = sys.maxsize
        for v in vli:
            d = self.dist(v.centroid())
            if d < mindist:
                mindist = d
        return mindist

    @lazy_property
    def dist_to_nearest_vesicle_border(self):
        """Return distance to nearest vesicle border.
           If self is within a vesicle, the distance is negative.
           If self is within a vesicle, the distance to the border of
           that vesicle is returned, even if self is actually closer to
           the border of another vesicle. NB: If vesicles are allowed
           to overlap, a point may be within more than one vesicle.
           The distance returned is then that to the nearest border
           of the containing vesicles."""
        vli = []
        if self.opt.prioritize_lumen:
            vli = [v for v in self.profile.vli if self.is_within_polygon(v)]
        if not vli:
            vli = self.profile.vli
        mindist = sys.maxsize
        minv = []
        for v in vli:
            d = self.perpend_dist_closed_path(v)
            if d < mindist:
                mindist = d
                minv = v
        if self.is_within_polygon(minv):
            mindist = -mindist
        return mindist

    def get_nearest_neighbour(self, pointli):
        """Determine distance to nearest neighbour."""
        # Assumes that only valid (projectable, within shell etc) points
        # are in pointli
        mindist = float(sys.maxsize)
        #minp = Point()
        for p in pointli:
            # I will instead exclude non-desired points from the
            # supplied point list *before* calling this function
            #if p is not self and p.isAssociatedWithProfile:
            if p is not self:
                d = self.dist(p)
                if d < mindist:
                    mindist = d
                    #minp = p
        if not mindist < float(sys.maxsize):
            return None
        else:
            self.nearest_neighbour_dist = mindist
            #self.nearest_neighbour_point = minp
            return self.nearest_neighbour_dist

    def get_nearest_lateral_neighbour(self, pointli):
        """Determine distance along profile border to nearest neighbour."""
        # Assumes that only valid (projectable, within shell etc) points
        # are in pointli
        mindist = float(sys.maxsize)
        minp = Point()
        for p in pointli:
            if p is not self:
                d = self.lateral_dist_to_point(p, self.profile.path)
                if d < mindist:
                    mindist = d
                    minp = p
        if not mindist < float(sys.maxsize):
            return None
        else:
            self.nearest_lateral_neighbour_dist = mindist
            self.nearest_lateral_neighbour_point = minp
            return self.nearest_lateral_neighbour_dist


class Vesicle(geometry.SegmentedPath):
    def __init__(self, pointlist=None, profile=None):
        if pointlist is None:
            pointlist = []
        self.profile = profile
        if self.profile is not None:
            self.opt = self.profile.opt
        else:
            self.opt = None
        geometry.SegmentedPath.__init__(self, pointlist)

    @lazy_property
    def center(self):
        return self.centroid()

    @lazy_property
    def circularity(self):
        """ Uses the same definition as ImageJ """
        return 4 * math.pi * (self.area() / math.pow(self.perimeter(), 2))

    @lazy_property
    def convexity(self):
        """ Return the ratio of the area of the vesicles and the area
            of its convex hull
        """
        return self.area() / self.convex_hull().area()

    @lazy_property
    def diameter_of_equal_area_circle(self):
        """ Return the diameter of a circle with the same area as self """
        return 2 * math.sqrt(self.area() / math.pi)

    @lazy_property
    def diameter_of_equal_perimeter_circle(self):
        """ Return the diameter of a circle with the same perimeter as self """
        return self.perimeter() / math.pi


class ProfileBorder(geometry.SegmentedPath):
    def __init__(self, pointlist=None, profile=None):
        if pointlist is None:
            pointlist = []
        geometry.SegmentedPath.__init__(self, pointlist)
        self.profile = profile


class PointList(list):
    def __init__(self, pointli, ptype, profile):
        super().__init__()
        try:
            self.extend([Point(p.x, p.y, ptype, profile) for p in pointli])
        except (AttributeError, IndexError):
            raise TypeError('not a list of Point elements')


class Profile:
    def __init__(self, inputfn, opt):
        self.inputfn = inputfn
        self.src_img = None
        self.outputfn = ""
        self.opt = opt
        self.pli = []
        self.vli = []
        self.holeli = []
        self.randomli = []
        self.interdistlis = {}
        self.interlatdistlis = {}
        self.n_discarded = {"particle": 0, "random": 0}
        self.comment = ""
        self.pixelwidth = None
        self.metric_unit = ""
        self.posloc = geometry.Point()
        self.negloc = geometry.Point()
        self.feret = None
        self.warnflag = False
        self.errflag = False

    def process(self):
        """ Parse profile data from a file and determine distances
        """
        try:
            self.__parse()
            self.__check_paths()
            sys.stdout.write("Processing profile...\n")
            self.__compute_stuff()
            # self.posloc = self.path.centroid()
            # if not self.posloc.is_within_polygon(self.path):
            #     self.negloc = self.posloc
            #     self.posloc = geometry.Point()
            if self.opt.determine_interparticle_dists:
                self.__determine_interdistlis("particle")
            if self.opt.determine_intervesicle_dists:
                self.__determine_interdistlis("vesicle")
            if self.opt.stop_requested:
                return
            sys.stdout.write("Done.\n")
        except ProfileError as err:
            sys.stdout.write("Error: %s\n" % err.msg)
            self.errflag = True

    @lazy_property
    def area(self):
        """Determine area of profile, excluding holes"""
        tot_hole_area = sum([h.area() for h in self.holeli])
        return self.path.area() - tot_hole_area

    def contains(self, p):
        """Determine if a point is inside profile, excluding holes."""
        if not p:
            return None
        return p.is_within_profile(self)

    def __compute_stuff(self):
        __ = self.area  # Force computation here
        self.perimeter = self.path.perimeter()
        self.feret = self.path.feret_diameter()
        for p in self.pli:
            p.determine_stuff()
        self.pli = [p for p in self.pli if not p.discard]
        for p in self.randomli:
            p.determine_stuff()
        self.randomli = [p for p in self.randomli if not p.discard]
        for ptype in ('particle', 'random'):
            if ptype == 'random' and not self.opt.use_random:
                continue
            ptypestr = 'particles' if ptype == 'particle' else ptype + ' points'
            sys.stdout.write("  Number of %s discarded: %d\n"
                             % (ptypestr, self.n_discarded[ptype]))

    def __determine_interdistlis(self, comp):
        if comp == 'particle':
            compli = self.pli
        elif comp == 'vesicle':
            compli = [Point(v.centroid(), profile=self) for v in self.vli]
        else:
            return
        rel_dict = self.opt.__dict__['inter%s_relations' % comp]
        dist_mode = self.opt.__dict__['inter%s_dist_mode' % comp]
        shortest_dist = self.opt.__dict__['inter%s_shortest_dist' % comp]
        lateral_dist = self.opt.__dict__['inter%s_lateral_dist' % comp]
        if True not in [val for key, val in rel_dict.items()]:
            return
        sys.stdout.write("Determining interparticle distances...\n")
        if rel_dict['%s - %s' % (comp, comp)]:
            self.interdistlis['%s - %s' % (comp, comp)] = \
                self.__get_same_interdistances(compli, dist_mode, shortest_dist, lateral_dist)
        if self.opt.use_random and rel_dict['random - %s' % comp]:
            self.interdistlis['random - %s' % comp] = \
                self.__get_interdistances(self.randomli, compli, dist_mode, shortest_dist,
                                          lateral_dist)
        if self.opt.use_random and rel_dict['%s - random' % comp]:
            self.interdistlis['%s - random' % comp] = \
                self.__get_interdistances(compli, self.randomli, dist_mode, shortest_dist,
                                          lateral_dist)

    def __get_same_interdistances(self, pointli, dist_mode, shortest_dist,
                                  lateral_dist):
        distdict = {'shortest': [], 'lateral': []}
        for i in range(0, len(pointli)):
            if self.opt.stop_requested:
                return [], []
            if dist_mode == 'all':
                for j in range(i + 1, len(pointli)):
                    if shortest_dist:
                        distdict['shortest'].append(pointli[i].dist(pointli[j]))
                    if lateral_dist:
                        distdict['lateral'].append(pointli[i].lateral_dist_to_point(pointli[j], 
                                                                                    self.path))
            elif dist_mode == 'nearest neighbour':
                if shortest_dist:
                    distdict['shortest'].append(pointli[i].get_nearest_neighbour(pointli))
                if lateral_dist:
                    distdict['lateral'].append(pointli[i].get_nearest_lateral_neighbour(pointli))
        distdict['shortest'] = [d for d in distdict['shortest'] if d is not None]
        distdict['lateral'] = [d for d in distdict['lateral'] if d is not None]
        return distdict

    def __get_interdistances(self, pointli, pointli2, dist_mode, shortest_dist, lateral_dist):
        if pointli2 is None:
            pointli2 = []
        distdict = {'shortest': [], 'lateral': []}
        for i, p in enumerate(pointli):
            if self.opt.stop_requested:
                return [], []
            if dist_mode == 'all':
                for p2 in pointli2:
                    if shortest_dist:
                        distdict['shortest'].append(p.dist(p2))
                    if lateral_dist:
                        distdict['lateral'].append(p.lateral_dist_to_point(p2, self.path))
            elif dist_mode == 'nearest neighbour':
                if shortest_dist:
                    distdict['shortest'].append(p.get_nearest_neighbour(pointli2))
                if lateral_dist:
                    distdict['lateral'].append(p.get_nearest_lateral_neighbour(pointli2))
        distdict['shortest'] = [d for d in distdict['shortest'] if d is not None]
        distdict['lateral'] = [d for d in distdict['lateral'] if d is not None]
        return distdict

    def __parse(self):
        """Parse profile data from input file."""
        sys.stdout.write("\nParsing '%s':\n" % self.inputfn)
        li = file_io.read_file(self.inputfn)
        if not li:
            raise ProfileError(self, "Could not open input file")
        while li:
            s = li.pop(0).replace('\n', '').strip()
            if s.split(' ')[0].upper() == "IMAGE":
                self.src_img = s.split(' ')[1]
            elif s.split(' ')[0].upper() == "PROFILE_ID":
                try:
                    self.id = s.split(' ')[1]
                except IndexError:
                    profile_warning(self, "Profile id not defined or invalid")
            elif s.split(' ')[0].upper() == "COMMENT":
                try:
                    self.comment = s.split(' ', 1)[1]
                except IndexError:
                    self.comment = ''
            elif s.split(' ')[0].upper() == "PIXELWIDTH":
                try:
                    self.pixelwidth = float(s.split(' ')[1])
                    self.metric_unit = s.split(' ')[2]
                except (IndexError, ValueError):
                    raise ProfileError(self, "PIXELWIDTH is not a valid number")
            elif s.upper() == "VESICLE":
                self.vli.append(Vesicle(self.__get_coords(li, 'vesicle'), self))
            elif s.upper() == "PROFILE_BORDER":
                self.path = ProfileBorder(self.__get_coords(li, 'path'), self)
            elif s.upper() in ("PROFILE_HOLE", "HOLE"):
                self.holeli.append(geometry.SegmentedPath(self.__get_coords(li, 'hole')))
            elif s.upper() in ("POINTS", "PARTICLES"):
                self.pli = PointList(self.__get_coords(li, 'particle'), 'particle', self)
            elif s.upper() == "RANDOM_POINTS":
                self.randomli = PointList(self.__get_coords(li, 'random'), 'random', self)
            elif s[0] != "#":  # unless specifically commented out
                profile_warning(self, "Unrecognized string '" + s + "' in input file")
        # Now, let's see if everything was found
        self.__check_parsed_data()

    def __check_parsed_data(self):
        """See if the profile data was parsed correctly, and print info
        on the parsed data to stdout.
        """
        self.__check_var_default('src_img', "Source image", "N/A")
        self.__check_var_default('id', "Profile id", "N/A")
        self.__check_var_default('comment', "Comment", "")
        self.__check_var_val('metric_unit', "Metric unit", 'metric_unit')
        self.__check_required_var('pixelwidth', "Pixel width", self.metric_unit)
        self.__check_list_var('path', 'Profile border', 'nodes', 2)
        self.__check_list_var('pli', 'Particles', '', 0)
        self.__check_table_var('vli', "Vesicle", "Vesicles", 1, 2)
        self.__check_table_var('holeli', "Hole", "Holes", 0, 2)
        self.__check_var_exists('randomli', "Random points", 'use_random')

    def __check_required_var(self, var_to_check, var_str, post_str):
        """Confirm that self has a required variable; else, raise
        ProfileError.
        """
        if not self.__dict__[var_to_check]:
            raise ProfileError(self, "%s not found in input file" % var_str)
        else:
            sys.stdout.write("  %s: %s %s\n" % (var_str, self.__dict__[var_to_check], post_str))

    @staticmethod
    def __check_list_len(var, min_len):
        """Return True if var is a list and has at least min_len
        elements, else False.
        """
        return isinstance(var, list) and len(var) >= min_len

    def __check_list_var(self, var_to_check, var_str, post_str, min_len):
        """Confirms that self has a var_to_check that is a list and
        has at least min_len elements; if var_to_check does not exist
        and min_len <= 0, assigns an empty list to var_to_check. Else,
        raise a ProfileError.
        """
        if not self.__dict__[var_to_check]:
            if min_len > 0:
                raise ProfileError(self, "%s not found in input file" % var_str)
            else:
                self.__dict__[var_to_check] = []
        elif not self.__check_list_len(self.__dict__[var_to_check], min_len):
            raise ProfileError(self, "%s has too few coordinates" % var_str)
        if post_str != '':
            post_str = " " + post_str
        sys.stdout.write("  %s%s: %d\n" % (var_str, post_str, len(self.__dict__[var_to_check])))

    def __check_table_var(self, var_to_check, var_str_singular,
                          var_str_plural, min_len_1, min_len_2):
        """Confirms that var_to_check exists, is a list and has at
        least min_len_1 elements, and that each of these has at least
        min_len_2 subelements; if var_to_check does not exist and
        min_len_1 <= 0, assigns an empty list to var_to_check. Else,
        raise ProfileError.
        """
        if not self.__dict__[var_to_check]:
            if min_len_1 > 0:
                raise ProfileError(self, "%s not found in input file" % var_str_plural)
            else:
                self.__dict__[var_to_check] = []
        elif not self.__check_list_len(self.__dict__[var_to_check], min_len_1):
            raise ProfileError(self, "Too few %s found in input file" % var_str_plural.lower())
        else:
            for element in self.__dict__[var_to_check]:
                if not self.__check_list_len(element, min_len_2):
                    raise ProfileError(self, "%s has too few coordinates" % var_str_singular)
        sys.stdout.write("  %s: %d\n" % (var_str_plural, len(self.__dict__[var_to_check])))

    def __check_var_default(self, var_to_check, var_str, default=""):
        """Checks if var_to_check exists; if not, assign the default
        value to var_to_check. Never raises a ProfileError.
        """
        if not self.__dict__[var_to_check]:
            self.__dict__[var_to_check] = default
        sys.stdout.write("  %s: %s\n" % (var_str, self.__dict__[var_to_check]))

    def __check_var_exists(self, var_to_check, var_str, optflag):
        """Checks for consistency between profiles with respect to the
        existence of var_to_check (i.e., var_to_check must be present
        either in all profiles or in none).

        If optflag is not set (i.e., this is the first profile), then
        set optflag to True or False depending on the existence of
        var_to_check. If optflag is already set (for consequent
        profiles), var_to_check must (if optflag is True) or must not
        (if optflag is False) exist. If not so, raise ProfileError.
        """
        if not hasattr(self.opt, optflag):
            if self.__dict__[var_to_check]:
                self.opt.__dict__[optflag] = True
            else:
                self.opt.__dict__[optflag] = False
        if self.opt.__dict__[optflag]:
            if self.__dict__[var_to_check]:
                sys.stdout.write("  %s: yes\n" % var_str)
            else:
                raise ProfileError(self, "%s expected but not found in input file" % var_str)
        elif self.__dict__[var_to_check]:
            raise ProfileError(self, "%s found but not expected" % var_str)
        else:
            sys.stdout.write("  %s: no\n" % var_str)

    def __check_var_val(self, var_to_check, var_str, optvar):
        """Checks for consistency between profiles with respect to the
        value of var_to_check (i.e., var_to_check must be present and
        have equal value in all profiles).

        If optvar is not set (i.e., this is the first profile), then
        set optflag to the value of var_to_check. If optvar is already
        set (for consequent profiles), the value of var_to_check must
        be equal to that of optvar. If not so, raise ProfileError.
        """
        if not self.__dict__[var_to_check]:
            raise ProfileError(self, "%s not found in input file" % var_str)
        if not hasattr(self.opt, optvar):
            self.opt.__dict__[optvar] = self.__dict__[var_to_check]
        elif self.__dict__[var_to_check] == self.opt.__dict__[optvar]:
            pass  # really no point in pointing out that it's ok
            # sys.stdout.write("  %s: %s\n"
            #                  % (var_str, parent.__dict__[var_to_check]))
        else:
            raise ProfileError(self, "%s value '%s'  differs from the value "
                                     "specified ('%s') in the first input file"
                               % (var_str, self.__dict__[var_to_check],
                                  self.opt.__dict__[optvar]))

    def __check_paths(self):
        """Check that the paths are valid, i e that:
           - profile border, vesicles and holes do not intersect with
             themselves;
           - vesicles and holes do not intersect with each other or
             with the profile border;
           - vesicles and holes are not within another vesicle or
             hole;
           - vesicles and holes are completely within the profile"""

        def check_simple(_path, s):
            if not _path.is_simple_polygon():
                raise ProfileError(self, "%s is invalid (crosses itself)" % s.capitalize())

        def check_in_profile(_path, s):
            if not _path.is_within_polygon(self.path):
                raise ProfileError(self, "%s is not within the profile " % s.capitalize())

        def check_overlap_own_type(pathli, s):
            for n, p in enumerate(pathli):
                for n2, p2 in enumerate(pathli):
                    if p is not p2 and p.overlaps_polygon(p2):
                        raise ProfileError(self, "%s %d overlaps with %s %d"
                                                 % (s.capitalize(), n + 1, s, n + n2 + 2))

        def check_overlap_other(pathli, otherpath, s, s_other):
            for n, p in enumerate(pathli):
                if p.crosses_polygon(otherpath):
                    raise ProfileError(self, "%s %d overlaps with %s"
                                             % (s.capitalize(), n + 1, s_other))
                if s_other != "profile border":
                    if p.is_within_polygon(otherpath):
                        raise ProfileError(self, "%s %d is within %s"
                                                 % (s.capitalize(), n + 1, s_other))
                    if otherpath.is_within_polygon(p):
                        raise ProfileError(self, "A %s is within %s %d" % (s, s_other, n + 1))

        check_simple(self.path, "profile border")
        check_overlap_other(self.vli, self.path, "vesicle", "profile border")
        check_overlap_other(self.holeli, self.path, "hole",
                            "profile border")
        for path in self.vli:
            check_simple(path, "vesicle")
            check_in_profile(path, "vesicle")
        for path in self.holeli:
            check_simple(path, "hole")
            check_in_profile(path, "hole")
        if not self.opt.allow_vesicle_overlap:
            check_overlap_own_type(self.vli, "vesicle")
        check_overlap_own_type(self.holeli, "hole")
        for hn, hole in enumerate(self.holeli):
            check_overlap_other(self.vli, hole, "vesicle", "hole %s" % (hn + 1))
        sys.stdout.write("  Paths are ok.\n")

    def __get_coords(self, strli, coord_type=""):
        """Pop point coordinates from list strli.

        When an element of strli is not a valid point, a warning is
        issued.
        """
        try:
            s = strli.pop(0).replace('\n', '').replace(' ', '').strip()
        except IndexError:
            return []
        pointli = []
        while s != 'END':
            try:
                p = geometry.Point(float(s.split(',')[0]), float(s.split(',')[1]))
                if pointli and (p == pointli[-1] or (coord_type == 'particle' and p in pointli)):
                    sys.stdout.write("Duplicate %s coordinates %s: skipping "
                                     "2nd instance\n" % (coord_type, p))
                else:
                    pointli.append(Point(p.x, p.y, ptype=coord_type))
            except ValueError:
                if s[0] != '#':
                    profile_warning(self, "'%s' not valid %s coordinates" % (s, coord_type))
                else:
                    pass
            try:
                s = strli.pop(0).replace('\n', '').strip()
            except IndexError:
                break
        # For some reason, sometimes the endnodes have the same coordinates;
        # in that case, delete the last endnode to avoid division by zero
        if (len(pointli) > 1) and (pointli[0] == pointli[-1]):
            del pointli[-1]
        return pointli
# end of class Profile


class OptionData:
    def __init__(self):
        self.input_file_list = []
        self.spatial_resolution = 25
        self.shell_width = 0
        self.outputs = {'profile summary': True,
                        'vesicle summary': True,
                        'particle summary': True,
                        'random summary': True,
                        'session summary': True}
        self.output_file_format = 'excel'
        self.output_filename_ext = '.xlsx'
        self.input_filename_ext = '.ves'
        self.output_filename_suffix = ''
        self.output_filename_other_suffix = ''
        self.output_filename_date_suffix = True
        self.csv_delimiter = 'comma'
        self.action_if_output_file_exists = 'overwrite'
        self.output_dir = ''
        self.allow_vesicle_overlap = True
        self.prioritize_lumen = False
        self.determine_particle_vesicle_center_dist = False
        self.determine_particle_vesicle_border_dist = False
        self.determine_intervesicle_dists = False
        self.intervesicle_dist_mode = 'nearest neighbour'
        self.intervesicle_relations = {'vesicle - vesicle': True,
                                       'vesicle - random': True,
                                       'random - vesicle': True}
        self.intervesicle_shortest_dist = True
        self.intervesicle_lateral_dist = True
        self.determine_interparticle_dists = False
        self.interparticle_dist_mode = 'nearest neighbour'
        self.interparticle_relations = {'particle - particle': True,
                                        'particle - random': True,
                                        'random - particle': True}
        self.interparticle_shortest_dist = True
        self.interparticle_lateral_dist = False
        self.metric_unit = None
        self.use_random = None
        self.particles_found = False
        self.stop_requested = False

    def reset(self):
        """ Resets all options to default, and removes those that are not
            set in __init__().
        """
        self.__dict__ = {}
        self.__init__()
# end of class OptionData


class ProfileError(Exception):
    def __init__(self, profile, msg):
        self.profile = profile
        self.msg = msg + "."


def profile_warning(profile, msg):
    """ Issue a warning
    """
    sys.stdout.write("Warning: %s.\n" % msg)
    profile.warnflag = True


def profile_message(msg):
    """ Show a message
    """
    sys.stdout.write("%s.\n" % msg)
