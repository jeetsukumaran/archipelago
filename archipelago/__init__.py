#! /usr/bin/env python

##############################################################################
## Copyright (c) 2014 Jeet Sukumaran.
## All rights reserved.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.
##     * Redistributions in binary form must reproduce the above copyright
##       notice, this list of conditions and the following disclaimer in the
##       documentation and/or other materials provided with the distribution.
##     * The names of its contributors may not be used to endorse or promote
##       products derived from this software without specific prior written
##       permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
## IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
## THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
## PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL JEET SUKUMARAN OR MARK T. HOLDER
## BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
## CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
## SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
## INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
## CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
## ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
## POSSIBILITY OF SUCH DAMAGE.
##
##############################################################################

import os

__project__ = "Archipelago"
__version__ = "0.1.0"
__archipelago_revision__ = None
__archipelago_description__ = None

ARCHIPELAGO_HOME = os.path.dirname(os.path.abspath(__file__))

def revision():
    global __archipelago_revision__
    if __archipelago_revision__ is None:
        from dendropy.utility import vcsinfo
        try:
            try:
                __homedir__ = os.path.dirname(os.path.abspath(__file__))
            except IndexError:
                __homedir__ = os.path.dirname(os.path.abspath(__file__))
        except OSError:
            __homedir__ = None
        except:
            __homedir__ = None
        __archipelago_revision__ = vcsinfo.Revision(repo_path=__homedir__)
    return __archipelago_revision__

def description():
    global __archipelago_description__
    if __archipelago_description__ is None:
        archipelago_revision = revision()
        if archipelago_revision.is_available:
            revision_text = " ({})".format(archipelago_revision)
        else:
            revision_text = ""
        __archipelago_description__  = "{} {}{}".format(__project__, __version__, revision_text)
    return __archipelago_description__
