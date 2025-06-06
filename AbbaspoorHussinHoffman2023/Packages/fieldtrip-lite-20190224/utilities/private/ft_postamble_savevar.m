% FT_POSTAMBLE_SAVEVAR is a helper script that optionally saves the output
% FieldTrip data structures to a *.mat file on disk. This is useful for
% batching and for distributed processing. This makes use of the
% cfg.outputfile variable.
%
% Use as
%   ft_postamble savevar data
%   ft_postamble savevar source mri
%
% See also FT_PREAMBLE, FT_POSTAMBLE, FT_POSTAMBLE_SAVEFIG

% Copyright (C) 2011-2019, Robert Oostenveld, DCCN
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% the output data should be saved to a MATLAB file
if (isfield(cfg, 'outputfile') && ~isempty(cfg.outputfile)) || exist('Fief7bee_reproducescript', 'var')
  
  if isfield(cfg, 'outputlock') && ~isempty(cfg.outputlock)
    mutexlock(cfg.outputlock);
  end
  
  if exist('Fief7bee_reproducescript', 'var')
    % write the function output variable(s) to a MATLAB file
    iW1aenge_now = datestr(now, 30);
    cfg.outputfile = {};
    for i=1:numel(iW1aenge_postamble)
      cfg.outputfile{i} = fullfile(Fief7bee_reproducescript, sprintf('%s_output_%s', iW1aenge_now, iW1aenge_postamble{i}));
    end
    % write a snippet of MATLAB code with the user-specified configuration and function call
    reproducescript(fullfile(Fief7bee_reproducescript, 'script.m'), cfg, isempty(iW1aenge_postamble))
  elseif (isfield(cfg, 'outputfile') && ~isempty(cfg.outputfile))
    % keep the output file as it is
  else
    % don't write to an output file
    cfg.outputfile = {};
  end
  
  if isequal(iW1aenge_postamble, {'varargout'}) && ~iscell(cfg.outputfile)
    % this should be a cell-array, oterwise it cannot be matched with varargout
    cfg.outputfile = {cfg.outputfile};
  end
  
  if iscell(cfg.outputfile)
    % iW1aenge_postamble is a cell-array containing the variable names
    if isequal(iW1aenge_postamble, {'varargout'})
      % the output is in varargout
      for tmpindx=1:length(cfg.outputfile)
        savevar(cfg.outputfile{tmpindx}, 'data', varargout{tmpindx});
      end % for
      clear tmpindx
    else
      % the output is in explicitly named variables
      for tmpindx=1:length(cfg.outputfile)
        savevar(cfg.outputfile{tmpindx}, iW1aenge_postamble{tmpindx}, eval(iW1aenge_postamble{tmpindx}));
      end % for
      clear tmpindx
    end
  else
    % iW1aenge_postamble{1} contains the name of the only variable
    savevar(cfg.outputfile, iW1aenge_postamble{1}, eval(iW1aenge_postamble{1}));
  end
  
  if isfield(cfg, 'outputlock') && ~isempty(cfg.outputlock)
    mutexunlock(cfg.outputlock);
  end
  
  if ~ft_nargout
    % do not return the output variable "ans"
    clear(iW1aenge_postamble{1});
  end
end

