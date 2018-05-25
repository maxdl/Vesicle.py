import configparser
import os
import os.path
import queue
import sys
import threading
import time
import traceback
import wx
from . import core
from . import file_io
from . import gui
from . import main
from . import stringconv
from . import version


class Frame(gui.MainFrame):
    def __init__(self, parent):
        gui.MainFrame.__init__(self, parent)
        self.SetTitle(version.title)
        self.SetIcon(wx.Icon(version.icon, wx.BITMAP_TYPE_ICO))
        self.set_win7_taskbar_icon()
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.set_input_file_list_ctrl_columns(self.InputFileListCtrl)
        dt = FileDropTarget(self)
        self.InputFileListCtrl.SetDropTarget(dt)
        self.opt = core.OptionData()
        self.configfn = os.path.normpath(os.path.expanduser('~/.%s.cfg' % version.title.lower()))
        self.log = None
        self.exitcode = None
        self.get_input_dir_from_config()
        self.load_options_from_config()
        self.set_options_in_ui()
        self.Fit()

    @staticmethod
    def set_win7_taskbar_icon():
        """ A hack to make the icon visible in the taskbar in Windows 7.
            From http://stackoverflow.com/a/1552105/674475.
        """
        if sys.platform == "win32":
            import ctypes
            appid = 'Vesicle'
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(appid)

    def OnAddFile(self, event):
        dlg = wx.FileDialog(self, "Choose a file", ".", "", "*%s"
                            % self.opt.input_filename_ext, wx.FD_MULTIPLE | wx.FD_CHANGE_DIR)
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.add_files(dlg.GetPaths())
        finally:
            dlg.Destroy()
    
    def OnRemoveFile(self, event):
        while 1:
            i = self.InputFileListCtrl.GetNextItem(-1, state=wx.LIST_STATE_SELECTED)
            if i == -1:
                break
            else:
                self.InputFileListCtrl.DeleteItem(i)
                
    def OnViewFile(self, event):
        if self.InputFileListCtrl.GetSelectedItemCount() == 0:
            self.show_warning("No file selected.")
            return
        elif self.InputFileListCtrl.GetSelectedItemCount() > 1:
            self.show_warning("You can only view one file at a time.")
            return
        i = self.InputFileListCtrl.GetNextItem(-1, state=wx.LIST_STATE_SELECTED)
        try:
            fn = os.path.join(self.InputFileListCtrl.GetItemText(i, 1),
                              self.InputFileListCtrl.GetItemText(i, 0))
        except IOError:
            self.show_error("Could not open file.")
            return
        dlg = ViewFileDialog(self, fn)
        try:
            dlg.ShowModal()
        finally:
            dlg.Destroy()

    def OnIntervesicleCheckbox(self, event):
        self.IntervesicleModeChoice.Enable(self.IntervesicleCheckBox.GetValue())
        self.IntervesicleModeLabel.Enable(self.IntervesicleCheckBox.GetValue())
        self.IntervesicleRelationsCheckListBox.Enable(self.IntervesicleCheckBox.GetValue())
        self.IntervesicleRelationsLabel.Enable(self.IntervesicleCheckBox.GetValue())
        self.IntervesicleShortestDistCheckBox.Enable(self.IntervesicleCheckBox.GetValue())
        self.IntervesicleLateralDistCheckBox.Enable(self.IntervesicleCheckBox.GetValue())

    def OnInterparticleCheckbox(self, event):
        self.InterparticleModeChoice.Enable(self.InterparticleCheckBox.GetValue())
        self.InterparticleModeLabel.Enable(self.InterparticleCheckBox.GetValue())
        self.InterparticleRelationsCheckListBox.Enable(self.InterparticleCheckBox.GetValue())
        self.InterparticleRelationsLabel.Enable(self.InterparticleCheckBox.GetValue())
        self.InterparticleShortestDistCheckBox.Enable(self.InterparticleCheckBox.GetValue())
        self.InterparticleLateralDistCheckBox.Enable(self.InterparticleCheckBox.GetValue())

    def OnOtherSuffixCheckBox(self, event):
        self.OtherSuffixTextCtrl.Enable(self.OtherSuffixCheckBox.GetValue())
        
    def OnSaveLogCheckBox(self, event):
        self.LogFilePickerCtrl.Enable(self.SaveLogCheckBox.GetValue())
        self.IfLogExistsRadioBox.Enable(self.SaveLogCheckBox.GetValue())
        
    def OnSetOptionsAsDefault(self, event):
        if self.save_options_to_config():
            self.StatusBar.SetStatusText("Current options saved to '%s'." % self.configfn)

    def OnStart(self, event):
        if self.InputFileListCtrl.GetItemCount() == 0:
            self.show_warning("No files to process.")
            return
        self.set_options_from_ui()
        if not self.set_log():
            return
        self.StatusBar.SetStatusText("Processing...")
        self.exitcode = 1
        event_type = ""
        msg = "Processing %s \n(File %d of %d)" % (
              os.path.basename(self.opt.input_file_list[0]), 1, len(self.opt.input_file_list))
        i = 0
        dlg = wx.ProgressDialog(version.title, msg,
                                len(self.opt.input_file_list) + 2,
                                parent=self,
                                style=wx.PD_ELAPSED_TIME | wx.PD_REMAINING_TIME | wx.PD_CAN_ABORT)
        pthread = ProcessThread(self.opt)
        pthread.start()
        while pthread.isAlive() or not pthread.process_queue.empty():
            if not pthread.process_queue.empty():
                (event_type, data) = pthread.process_queue.get()
                if event_type == "new_file":
                    i += 1
                    msg = "Processing %s \n(File %d of %d)" \
                        % (os.path.basename(data), i,
                           len(self.opt.input_file_list))
                if event_type == "saving_summaries":
                    i += 1
                    msg = "Saving summaries..."
                if event_type == "done":
                    i = len(self.opt.input_file_list) + 2
                    msg = "Done."
            self.log.update()
            if event_type == "done" and self.log.fn != "":
                self.StatusBar.SetStatusText("Logged to '" + self.log.fn + "'.")
            if not dlg.Update(i, msg)[0] and not self.opt.stop_requested:
                if self.yes_no_dialog("Abort process?"):
                    pthread.stop()
                    dlg.Hide()
                else:
                    dlg.Resume()
            if dlg.GetSize().GetWidth() < dlg.GetBestSize().GetWidth():
                dlg.SetSize((dlg.GetBestSize().GetWidth() + 20, dlg.GetBestSize().GetHeight()))
        if not pthread.error_queue.empty():
            exc_str = pthread.error_queue.get()
            if self.log.fn != "":
                self.StatusBar.SetStatusText("Logged to '" + self.log.fn + "'.")
            sys.stdout.write("\n*** %s session was unexpectedly aborted"
                             " at %s (local time). \n\nDetails:\n%s"
                             % (version.title, time.ctime(), exc_str))
            self.log.update()
            self.show_error("An unexpected error occurred while executing "
                            "%s - session aborted.\n\nDetails (also "
                            "sent to log):\n\n %s" % (version.title, exc_str))
            dlg.Destroy()
            return
        # Processing finished.
        self.log.update()
        if self.log.fn != "":
            self.StatusBar.SetStatusText("Logged to '" + self.log.fn + "'.")
        dlg.Destroy()
        if pthread.exitcode == 0:
            self.show_error("One or more errors occurred during processing. "
                            "See log for details.")
        elif pthread.exitcode == 2:
            self.show_warning("One or more warnings occurred during "
                              "processing. See log for details.")
        elif pthread.exitcode == 3:
            self.show_warning("Session aborted by user.")
            
    def OnAbout(self, event):
        dlg = AboutDialog(self)
        try:
            dlg.ShowModal()
        finally:
            dlg.Destroy()

    def OnClose(self, event):
        self.save_input_dir_to_config()
        sys.stdout = sys.__stdout__
        self.Destroy()

#
#   utilities
#
    def save_input_dir_to_config(self):
        config = configparser.ConfigParser()
        try:
            config.read(self.configfn)
        except (configparser.ParsingError, configparser.MissingSectionHeaderError):
            pass  # Silently suppress parsing errors at this stage
        if 'Previous session' not in config.sections():
            config['Previous session'] = {}
        config['Previous session']['input_dir'] = os.getcwd()
        try:
            with open(self.configfn, 'w') as f:
                config.write(f)
        except IOError:
            self.show_warning("Configuration file\n(%s)\ncould not be saved."
                              % self.configfn)

    def get_input_dir_from_config(self):
        config = configparser.ConfigParser()
        if not os.path.exists(self.configfn):
            return
        try:
            config.read(self.configfn)
        except (configparser.ParsingError, configparser.MissingSectionHeaderError):
            pass  # Silently suppress parsing errors at this stage
        try:
            inputdir = config['Previous session']['input_dir']
        except (configparser.NoSectionError, configparser.NoOptionError, KeyError):
            self.show_warning("Configuration file '%s' is invalid.\n Using "
                              "current working directory." % self.configfn)
            return
        try:
            if not os.path.isdir(inputdir):
                raise IOError
            os.chdir(inputdir)
        except (IOError, TypeError):
            self.show_warning("Invalid input directory' %s' in configuration "
                              "file '%s'.\n Using current working directory."
                              % (inputdir, self.configfn))

    def save_options_to_config(self):

        def set_option(option):
            config.set('Options', option, str(getattr(self.opt, option)))

        def set_dict_option(option):
            optdict = getattr(self.opt, option)
            for key, val in optdict.items():
                optstr = '.'.join([option, key.replace(' ', '_')])
                config.set('Options', optstr, str(val))

        self.set_options_from_ui()
        config = configparser.ConfigParser()
        try:
            config.read(self.configfn)
        except (configparser.ParsingError, configparser.MissingSectionHeaderError):
            pass  # Silently suppress parsing errors at this stage
        if not config.has_section('Options'):
            config.add_section('Options')
        set_option('output_file_format')
        set_option('csv_delimiter')
        set_option('action_if_output_file_exists')
        set_option('output_filename_date_suffix')
        set_option('spatial_resolution')
        set_option('shell_width')
        set_option('allow_vesicle_overlap')
        set_option('prioritize_lumen')
        set_option('determine_intervesicle_dists')
        set_option('determine_particle_vesicle_center_dists')
        set_option('determine_particle_vesicle_border_dists')
        set_option('intervesicle_dist_mode')
        set_option('intervesicle_shortest_dist')
        set_option('intervesicle_lateral_dist')
        set_dict_option('intervesicle_relations')
        set_option('determine_interparticle_dists')
        set_option('interparticle_dist_mode')
        set_option('interparticle_shortest_dist')
        set_option('interparticle_lateral_dist')
        set_dict_option('interparticle_relations')
        set_dict_option('outputs')
        try:
            with open(self.configfn, 'wb') as f:
                config.write(f)
        except IOError:
            self.show_warning("Configuration file\n(%s)\ncould not be saved." % self.configfn)

    def load_options_from_config(self):

        def show_invalid_option_warning(invalid_opt):
            self.show_warning("Invalid value '%s' for option '%s' in "
                              "configuration file '%s'.\nUsing default value."
                              % (getattr(self.opt, invalid_opt), invalid_opt,
                                 self.configfn))

        def check_str_option(opt, valid_strings=()):
            if getattr(self.opt, opt) not in valid_strings:
                show_invalid_option_warning(opt)
                setattr(self.opt, opt, getattr(defaults, opt))

        def check_int_option(opt, lower=None, upper=None):
            try:
                setattr(self.opt, opt,
                        stringconv.str_to_int(getattr(self.opt, opt), lower, upper))
            except ValueError:
                show_invalid_option_warning(opt)
                setattr(self.opt, opt, getattr(defaults, opt))

        def check_bool_option(opt):
            try:
                setattr(self.opt, opt, stringconv.str_to_bool(getattr(self.opt, opt)))
            except ValueError:
                show_invalid_option_warning(opt)
                setattr(self.opt, opt, getattr(defaults, opt))

        def check_bool_dict_option(opt):
            optdict = getattr(self.opt, opt)
            defaultdict = getattr(defaults, opt)
            for key, val in list(optdict.items()):
                optstr = '.'.join([opt, key.replace(" ", "_")])
                if key not in list(getattr(defaults, opt).keys()):
                    self.show_warning("Invalid option '%s' in configuration file '%s'."
                                      % (optstr, self.configfn))
                    del optdict[key]
                try:
                    optdict[key] = stringconv.str_to_bool(val)
                except ValueError:
                    self.show_warning("Invalid value '%s' for option '%s' in "
                                      "configuration file '%s'.\nUsing default "
                                      "value." % (val, optstr, self.configfn))
                    optdict[key] = defaultdict[key]


        config = configparser.ConfigParser()
        if not os.path.exists(self.configfn):
            return
        try:
            config.read(self.configfn)
        except (configparser.ParsingError, configparser.MissingSectionHeaderError):
            return     # Silently suppress parsing errors at this stage
        if 'Options' not in config.sections():
            return     # No options present in config file; silently use defaults
        defaults = core.OptionData()
        for option in config.options('Options'):
            if '.' in option:
                option_dict, option_key = option.split('.', 1)
                option_key = option_key.replace("_", " ")
                try:
                    getattr(self.opt,
                            option_dict)[option_key] = config.get('Options', option)
                except AttributeError:
                    pass   # So, attribute is invalid, but continue silently
            else:
                setattr(self.opt, option, config.get('Options', option))
        check_str_option('output_file_format', ('excel', 'csv'))
        check_str_option('csv_delimiter', ('comma', 'tab'))
        check_str_option('action_if_output_file_exists', ('enumerate', 'overwrite'))
        check_bool_option('output_filename_date_suffix')
        check_int_option('spatial_resolution', lower=0, upper=1000)
        check_int_option('shell_width', lower=0, upper=10000)
        check_bool_option('allow_vesicle_overlap')
        check_bool_option('prioritize_lumen')
        check_bool_option('determine_particle_vesicle_center_dists')
        check_bool_option('determine_particle_vesicle_border_dists')
        check_bool_option('determine_intervesicle_dists')
        check_str_option('intervesicle_dist_mode', ('nearest neighbour', 'all'))
        check_bool_dict_option('intervesicle_relations')
        check_bool_option('determine_interparticle_dists')
        check_str_option('interparticle_dist_mode', ('nearest neighbour', 'all'))
        check_bool_option('interparticle_shortest_dist')
        check_bool_option('interparticle_lateral_dist')
        check_bool_dict_option('interparticle_relations')
        check_bool_dict_option('outputs')

    def set_options_in_ui(self):
        # Vesicle options
        self.VesicleOverlapCheckBox.SetValue(self.opt.allow_vesicle_overlap)
        self.PrioritizeLumenCheckBox.SetValue(self.opt.prioritize_lumen)
        self.IntervesicleCheckBox.SetValue(self.opt.determine_intervesicle_dists)
        self.IntervesicleModeChoice.SetItems(['Nearest neighbour', 'All'])
        self.IntervesicleModeChoice.SetStringSelection(self.opt.intervesicle_dist_mode)
        # In order to get it sorted the way I want, I have to set this
        # explicitly rather than fetch from the dict
        self.IntervesicleRelationsCheckListBox.SetItems(['Vesicle - vesicle',
                                                         'Vesicle - random',
                                                         'Random - vesicle'])
        self.IntervesicleRelationsCheckListBox.SetCheckedStrings(
            [key.capitalize() for key in self.opt.intervesicle_relations
                if self.opt.intervesicle_relations[key] is True])
        self.IntervesicleShortestDistCheckBox.SetValue(self.opt.intervesicle_shortest_dist)
        self.InterparticleLateralDistCheckBox.SetValue(self.opt.intervesicle_lateral_dist)
        self.IntervesicleModeChoice.Enable(self.IntervesicleCheckBox.GetValue())
        self.IntervesicleModeLabel.Enable(self.IntervesicleCheckBox.GetValue())
        self.IntervesicleRelationsCheckListBox.Enable(self.IntervesicleCheckBox.GetValue())
        self.IntervesicleRelationsLabel.Enable(self.IntervesicleCheckBox.GetValue())
        self.IntervesicleShortestDistCheckBox.Enable(self.IntervesicleCheckBox.GetValue())
        self.IntervesicleLateralDistCheckBox.Enable(self.IntervesicleCheckBox.GetValue())
        # Particle options
        self.SpatResSpinCtrl.SetValue(self.opt.spatial_resolution)
        self.ShellWidthSpinCtrl.SetValue(self.opt.shell_width)
        self.InterparticleCheckBox.SetValue(self.opt.determine_interparticle_dists)
        self.InterparticleModeChoice.SetItems(['Nearest neighbour', 'All'])
        self.InterparticleModeChoice.SetStringSelection(self.opt.interparticle_dist_mode)
        # In order to get it sorted the way I want, I have to set this
        # explicitly rather than fetch from the dict
        self.InterparticleRelationsCheckListBox.SetItems(['Particle - particle',
                                                       'Particle - random',
                                                       'Random - particle'])
        self.InterparticleRelationsCheckListBox.SetCheckedStrings(
            [key.capitalize() for key in self.opt.interparticle_relations
             if self.opt.interparticle_relations[key] is True])
        self.InterparticleShortestDistCheckBox.SetValue(self.opt.interparticle_shortest_dist)
        self.InterparticleLateralDistCheckBox.SetValue(self.opt.interparticle_lateral_dist)
        self.InterparticleModeChoice.Enable(self.InterparticleCheckBox.GetValue())
        self.InterparticleModeLabel.Enable(self.InterparticleCheckBox.GetValue())
        self.InterparticleRelationsCheckListBox.Enable(self.InterparticleCheckBox.GetValue())
        self.InterparticleRelationsLabel.Enable(self.InterparticleCheckBox.GetValue())
        self.InterparticleShortestDistCheckBox.Enable(self.InterparticleCheckBox.GetValue())
        self.InterparticleLateralDistCheckBox.Enable(self.InterparticleCheckBox.GetValue())
        # Output options
        self.OutputCheckListBox.SetCheckedStrings(
            [key.capitalize() for key in self.opt.outputs
                if self.opt.outputs[key] is True])
        if self.opt.output_file_format == 'excel':
            self.OutputFormatRadioBox.SetStringSelection('Excel')
        elif self.opt.csv_delimiter == 'comma':
            self.OutputFormatRadioBox.SetStringSelection('Comma-delimited text')
        else:
            self.OutputFormatRadioBox.SetSetStringSelection('Tab-delimited text')
        self.IfOutputExistsRadioBox.SetStringSelection(
            self.opt.action_if_output_file_exists.capitalize())
        self.DateSuffixCheckBox.SetValue(self.opt.output_filename_date_suffix)
        self.OtherSuffixCheckBox.SetValue(self.opt.output_filename_other_suffix != '')
        self.OtherSuffixTextCtrl.SetValue(self.opt.output_filename_other_suffix)
        self.OtherSuffixTextCtrl.Enable(self.OtherSuffixCheckBox.GetValue())
        # Log options
        self.LogFilePickerCtrl.SetPath(version.title + '.log')
 
    def set_options_from_ui(self):
        # Input files
        self.opt.input_file_list = []
        for n in range(0, self.InputFileListCtrl.GetItemCount()):
            self.opt.input_file_list.append(os.path.join(
                self.InputFileListCtrl.GetItemText(n, 1),
                self.InputFileListCtrl.GetItemText(n, 0)))
        # Vesicle options
        self.opt.allow_vesicle_overlap = self.VesicleOverlapCheckBox.GetValue()
        self.opt.prioritize_lumen = self.PrioritizeLumenCheckBox.GetValue()
        self.VesicleOverlapCheckBox.SetValue(self.opt.allow_vesicle_overlap)
        self.opt.determine_intervesicle_dists = self.IntervesicleCheckBox.GetValue()        
        for key in self.opt.intervesicle_relations:
            if key.capitalize() in self.IntervesicleRelationsCheckListBox.GetCheckedStrings():
                self.opt.intervesicle_relations[key] = True
            else:
                self.opt.intervesicle_relations[key] = False
        if True not in self.opt.intervesicle_relations.values():
            self.opt.determine_intervesicle_dists = False
        self.opt.intervesicle_dist_mode = self.IntervesicleModeChoice.GetStringSelection().lower()
        self.opt.intervesicle_shortest_dist = self.IntervesicleShortestDistCheckBox.GetValue()
        self.opt.intervesicle_lateral_dist = self.IntervesicleLateralDistCheckBox.GetValue()
        # Particle options
        self.opt.spatial_resolution = int(self.SpatResSpinCtrl.GetValue())
        self.opt.shell_width = int(self.ShellWidthSpinCtrl.GetValue())        
        self.opt.determine_interparticle_dists = self.InterparticleCheckBox.GetValue()
        for key in self.opt.interparticle_relations:
            if key.capitalize() in self.InterparticleRelationsCheckListBox.GetCheckedStrings():
                self.opt.interparticle_relations[key] = True
            else:
                self.opt.interparticle_relations[key] = False
        if True not in self.opt.interparticle_relations.values():
            self.opt.determine_interparticle_dists = False
        self.opt.interparticle_dist_mode = self.InterparticleModeChoice.GetStringSelection().lower()
        self.opt.interparticle_shortest_dist = self.InterparticleShortestDistCheckBox.GetValue()
        self.opt.interparticle_lateral_dist = self.InterparticleLateralDistCheckBox.GetValue()
        # Output options
        for key in self.opt.outputs:
            if key.capitalize() in self.OutputCheckListBox.GetCheckedStrings():
                self.opt.outputs[key] = True
            else:
                self.opt.outputs[key] = False
        if self.OutputFormatRadioBox.GetStringSelection() == "Excel":
            self.opt.output_file_format = 'excel'
            self.opt.output_filename_ext = '.xlsx'
        elif self.OutputFormatRadioBox.GetStringSelection() == "Comma-delimited text":
            self.opt.output_file_format = 'csv'
            self.opt.output_filename_ext = '.csv'
            self.opt.csv_delimiter = "comma"
        elif self.OutputFormatRadioBox.GetStringSelection() == "Tab-delimited text":
            self.opt.output_file_format = 'csv'
            self.opt.output_filename_ext = '.csv'
            self.opt.csv_delimiter = "tab"
        self.opt.action_if_output_file_exists = \
            self.IfOutputExistsRadioBox.GetStringSelection().lower()
        self.opt.output_filename_date_suffix = self.DateSuffixCheckBox.GetValue()
        if self.OtherSuffixCheckBox.GetValue():
            self.opt.output_filename_other_suffix = self.OtherSuffixTextCtrl.GetValue()
        self.get_output_dir()

    def get_input_dir(self):
        for f in self.opt.input_file_list:
            if os.path.dirname(f):
                return os.path.dirname(f)
        return ""

    def get_output_dir(self):
        self.opt.output_dir = os.path.join(self.get_input_dir() or os.getcwd(), "out")
        if not os.path.isdir(self.opt.output_dir):
            os.mkdir(self.opt.output_dir)

    def add_files(self, fli):
        if len(fli) == 0:
            return
        c = self.InputFileListCtrl.GetItemCount()
        n = 0
        fn = ""
        for fn in fli:
            if (os.path.isfile(fn) and
                    os.path.splitext(fn)[1] == self.opt.input_filename_ext):
                self.InputFileListCtrl.InsertItem(c + n, os.path.basename(fn))
                self.InputFileListCtrl.SetItem(c + n, 1, os.path.dirname(fn))
                n += 1
            elif os.path.isdir(fn):
                for fn2 in os.listdir(fn):
                    if (os.path.isfile(os.path.join(fn, fn2)) and
                            os.path.splitext(fn2)[1] ==
                            self.opt.input_filename_ext):
                        self.InputFileListCtrl.InsertItem(c + n, fn2)
                        self.InputFileListCtrl.SetItem(c + n, 1, fn)
                        n += 1
        if n > 0:
            self.InputFileListCtrl.SetColumnWidth(0, -1)
            self.InputFileListCtrl.SetColumnWidth(1, -1)
        elif os.path.isdir(fn):
            self.show_warning("No files with '%s' extension found in folder(s)."
                              % self.opt.input_filename_ext)
        else:
            self.show_warning("Input files must have a '%s' extension."
                              % self.opt.input_filename_ext)

    @staticmethod
    def set_input_file_list_ctrl_columns(parent):
        parent.InsertColumn(col=0, format=wx.LIST_FORMAT_LEFT, heading='Name', width=-1)
        parent.InsertColumn(col=1, format=wx.LIST_FORMAT_LEFT, heading='Path', width=-1)

    def set_log(self):  # hm can't I simplify this?
        if self.SaveLogCheckBox.GetValue():
            mode = self.IfLogExistsRadioBox.GetStringSelection()
            logfn = self.LogFilePickerCtrl.GetPath()
            if os.path.dirname(logfn) == "":
                logfn = os.path.join(self.opt.output_dir, logfn)
            try:
                if os.path.exists(logfn):
                    if self.IfLogExistsRadioBox.GetStringSelection() == "Enumerate":
                        logfn = file_io.enum_filename(logfn, 2)
                    else:
                        f = open(logfn, 'a')
                        f.close()
                # ok, so file doesn't exist but check if name is valid
                else:
                    f = open(logfn, 'w')
                    f.close()
            except IOError:
                self.show_error("Could not write to log file. Please choose another filename.")
                return 0
            self.log = LogQueue(self, logfn, self.LogTextCtrl, mode)
        else:
            self.log = LogQueue(self, "", self.LogTextCtrl, "")
        sys.stdout = self.log
        return 1

    def show_warning(self, s):
        dlg = wx.MessageDialog(self, s, version.title, wx.OK | wx.ICON_EXCLAMATION)
        try:
            dlg.ShowModal()
        finally:
            dlg.Destroy()

    def show_error(self, s):
        dlg = wx.MessageDialog(self, s, version.title, wx.OK | wx.ICON_HAND)
        try:
            dlg.ShowModal()
        finally:
            dlg.Destroy()

    def yes_no_dialog(self, s):
        dlg = wx.MessageDialog(self, s, version.title, wx.YES_NO | wx.ICON_QUESTION | wx.NO_DEFAULT)
        try:
            pressed = dlg.ShowModal()
        finally:
            dlg.Destroy()
        if pressed == wx.ID_YES:
            return True
        return False

    def yes_no_warn_dialog(self, s):
        dlg = wx.MessageDialog(self, s, version.title, wx.YES_NO | wx.ICON_WARNING | wx.NO_DEFAULT)
        try:
            pressed = dlg.ShowModal()
        finally:
            dlg.Destroy()
        if pressed == wx.ID_YES:
            return True
        return False


class ProcessThread(threading.Thread):
    def __init__(self, opt):
        threading.Thread.__init__(self)
        self.opt = opt
        self.process_queue = queue.Queue()
        self.error_queue = queue.Queue(1)
        self.opt.stop_requested = False
        self.exitcode = None

    def stop(self):
        self.opt.stop_requested = True

    # noinspection PyBroadException
    def run(self):
        try:
            self.exitcode = main.main_proc(self)
        except:  # yes, I do want to catch everything
            exc_str = "".join(traceback.format_exception(*sys.exc_info()))
            self.error_queue.put(exc_str)


class LogQueue:
    def __init__(self, parent, fn, win, mode):
        self.parent = parent
        self.fn = fn
        self.win = win
        self.q = queue.Queue()
        if self.fn != "":
            self.errstr = "* Error: could not write to log file: %s\n" % self.fn
            if mode == 'Append':
                try:
                    f = open(self.fn, "a")
                    f.close()
                except IOError:
                    try:
                        f = open(self.fn, "w")
                        f.close()
                    except IOError:
                        sys.stderr.write(self.errstr)
                        self.fn = ""
            elif mode == 'Overwrite' or mode == 'Enumerate':
                try:
                    f = open(self.fn, "w")
                    f.close()
                except IOError:
                    sys.stderr.write(self.errstr)
                    self.fn = ""

    def write(self, s):
        self.q.put(s)

    def update(self):
        while not self.q.empty():
            s = self.q.get()
            self.win.write(s)
            self.parent.Update()
            if self.fn != "":
                try:
                    f = open(self.fn, "a")
                    try:
                        f.write(s)
                    finally:
                        f.close()
                except IOError:
                    sys.stderr.write(self.errstr)


class FileDropTarget(wx.FileDropTarget):

    def __init__(self, parent):
        wx.FileDropTarget.__init__(self)
        self.parent = parent

    # noinspection PyMethodOverriding
    def OnDropFiles(self, x, y, fli):
        self.parent.add_files(fli)
        return True


class AboutDialog(gui.AboutDialog):

    def __init__(self, parent):
        gui.AboutDialog.__init__(self, parent)
        self.TitleLabel.SetLabel(version.title)
        self.IconBitmap.SetBitmap(wx.Bitmap(version.icon, wx.BITMAP_TYPE_ANY))
        self.VersionLabel.SetLabel("Version %s" % version.version)
        self.LastModLabel.SetLabel("Last modified %s %s, %s." % version.date)
        self.CopyrightLabel.SetLabel("Copyright 2001-%s %s." % (version.date[2], version.author))
        self.LicenseLabel.SetLabel("Released under the terms of the MIT"
                                   " license.")
        self.EmailHyperlink.SetLabel("%s" % version.email)
        self.EmailHyperlink.SetURL("mailto://%s" % version.email)
        self.WebHyperlink.SetLabel("http://%s" % version.homepage)
        self.WebHyperlink.SetURL("http://%s" % version.homepage)
        self.SetIcon(wx.Icon(version.icon, wx.BITMAP_TYPE_ICO))
        self.Fit()

    def OnClose(self, event):
        self.Destroy()


class ViewFileDialog(gui.ViewFileDialog):
    def __init__(self, parent, fn):
        gui.ViewFileDialog.__init__(self, parent)
        self.SetIcon(wx.Icon(version.icon, wx.BITMAP_TYPE_ICO))
        try:
            self.SetTitle(os.path.basename(fn))
            f = open(fn, "r")
            try:
                for s in f.readlines():
                    self.ViewFileTextCtrl.AppendText(s)
            finally:
                f.close()
        except IOError:
            parent.show_error("Could not open file.")
            self.Close()

    def OnClose(self, event):
        self.Destroy()
