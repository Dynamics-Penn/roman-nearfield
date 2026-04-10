"""
latex_outputs.py
----------
Functions for printing out tables for Latex. Claude helped with regex and debugging.
"""

import numpy as np
import astropy.coordinates
import astropy.table as table
from astropy import units as u
import re
from collections import Counter


def make_defcitealias(bibcodes, nicknames):
    if isinstance(bibcodes, np.ndarray):
        bibcodes = bibcodes.tolist()
    if isinstance(nicknames, np.ndarray):
        nicknames = nicknames.tolist()
    if len(bibcodes) != len(nicknames):
        raise ValueError("bibcodes and nicknames must be the same length")
    return [rf"\defcitealias{{{bib}}}{{{nick}}}" for bib, nick in zip(bibcodes, nicknames)]


def make_citepalias(bibcodes_to_cite, bibcodes, nicknames):
    if isinstance(bibcodes, np.ndarray):
        bibcodes = bibcodes.tolist()
    if isinstance(nicknames, np.ndarray):
        nicknames = nicknames.tolist()
    if isinstance(bibcodes_to_cite, np.ndarray):
        bibcodes_to_cite = bibcodes_to_cite.tolist()
    if len(bibcodes) != len(nicknames):
        raise ValueError("bibcodes and nicknames must be the same length")
    if isinstance(bibcodes_to_cite, str):
        bibcodes_to_cite = [bibcodes_to_cite]
    resolved = []
    for bibcode in bibcodes_to_cite:
        if bibcode not in bibcodes:
            raise ValueError(f"{bibcode!r} not found in bibcodes list")
        resolved.append(nicknames[bibcodes.index(bibcode)])
    return rf"\citepalias{{{', '.join(resolved)}}}"

def get_author(r):
    author = re.split(r'\d{4}', r)[0]     # everything before the year
    if author.lower().startswith('vande'):
        return 'vd' + author[5:]
    return author

def get_prefix(author, tl):
    if author.lower().startswith('vd'):
        return author[:tl+2]  # 'vd' counts as 2 extra characters
    return author[:tl]

def get_yearcode(r):
    match = re.search(r'\d{4}', r)
    if match is None:
        raise ValueError(f"No 4-digit year found in: {r!r}")
    return match.group()[2:]

def parse_references(mag_refs, dist_refs):

    mrl = mag_refs.filled().tolist()
    if np.any(mag_refs.mask):
        mrl.remove(mag_refs.fill_value)
    drl = dist_refs.filled().tolist()
    if np.any(dist_refs.mask):
        drl.remove(dist_refs.fill_value)
    allrefs = np.unique(mrl+drl)

    # pre-extract author strings and year codes separately

    authors = np.array([get_author(r) for r in allrefs])
    yearcodes = np.array([get_yearcode(r) for r in allrefs])

    # build reference codes
    tl = 1
    refcodes_list = [get_prefix(a, tl) + y for a, y in zip(authors, yearcodes)]

    # resolve duplicates
    un, num = np.unique(refcodes_list, return_counts=True)
    while un.shape[0] != len(refcodes_list):
        print('found duplicate codes:', un[num>1])
        counts = Counter(refcodes_list)
        dupes = [i for i, n in enumerate(refcodes_list) if counts[n] > 1]
        tl += 1
        for i in dupes:
            refcodes_list[i] = get_prefix(authors[i], tl) + yearcodes[i]
        un, num = np.unique(refcodes_list, return_counts=True)
    
    refcodes = np.array(refcodes_list)  # convert to array 

    #print out the preamble commands needed
    print('\n-----paste in preamble----\n')
    for cmd in make_defcitealias(allrefs, refcodes):
        print(cmd)
        
    return allrefs, refcodes

def print_latex_table_galaxies(dmin, dmax, combined_map=None, nbgs_all=None, nbg_coords_all=None):
    
    labels = np.array(["HLWAS WIDE", "HLWAS MEDIUM"])
    if combined_map is None:
        combined_map = get_footprint()
    if nbgs_all is None or nbg_coords_all is None:
        nbgs_all, nbg_coords_all = read_galaxies_all(dmax=1.05*dmax)
    dist = nbgs_all['distance']*u.kpc
    dsel = np.where((dist>dmin)*(dist<dmax))[0]
    idist = np.argsort(dist[dsel])

    mag_refs = np.unique(nbgs_all['ref_m_v'][dsel])
    dist_refs = np.unique(nbgs_all['ref_distance'][dsel])
    allrefs, refcodes = parse_references(mag_refs, dist_refs)

    print('\n-----paste in table----\n')
    for i in dsel[idist]:
        sat = nbg_coords_all[i]
        ipix = hp.ang2pix(hp.get_nside(combined_map),sat.galactic.l.value, sat.galactic.b.value,lonlat=True)
        if not combined_map[ipix]==hp.UNSEEN:
            which_s = labels[int(combined_map[ipix])-1][6:].lower()
            if sat.galactic.b.value > 0: 
                which_s += ' N'
            else:
                which_s += ' S'
            smass_oom = np.trunc(nbgs_all['mass_stellar'][i])
            smass_absc = 10**(nbgs_all['mass_stellar'][i] - smass_oom )
            refs = np.unique([nbgs_all['ref_m_v'][i],nbgs_all['ref_distance'][i]]).tolist()
            out = nbgs_all['name'][i]+' & '
            out += '${0:0.2f} \\times 10^{{{1}}}$'.format(smass_absc, int(smass_oom))
            out += ' & '
            out += '{0:=02.0f} {1:=02.0f} {2:=02.1f}'.format(nbg_coords_all[i].ra.hms.h,nbg_coords_all[i].ra.hms.m,nbg_coords_all[i].ra.hms.s) 
            out += ' & '
            out += '{0:+d} {1:=02.0f} {2:=02.1f}'.format(int(nbg_coords_all[i].dec.dms.d),np.abs(nbg_coords_all[i].dec.dms.m),np.abs(nbg_coords_all[i].dec.dms.s))
            out += ' & '
            # lat = nbg_coords_all[i].galactic.l.wrap_at(180*u.deg)
            # out += '{0:+03d} {1:=02.0f} {2:=02.1f}'.format(int(lat.dms.d),np.abs(lat.dms.m),np.abs(lat.dms.s))
            # out += ' & '
            # out += '{0:+03d} {1:=02.0f} {2:=02.1f}'.format(int(nbg_coords_all[i].galactic.b.dms.d),np.abs(nbg_coords_all[i].galactic.b.dms.m),np.abs(nbg_coords_all[i].galactic.b.dms.s))
            # out += ' & '
            out += '{0:0.1f}'.format(nbgs_all['M_V'][i])
            out += ' & '
            out += '{0:0.3f}'.format((dist[i].to(u.Mpc)).value)
            out += ' & '
            out += '{0:0.1f}'.format(nbgs_all['apparent_magnitude_v'][i])
            out += ' & ' + which_s
            out += ' & ' + make_citepalias(refs,allrefs,refcodes)
            out += ' \\\\'

            print(out)