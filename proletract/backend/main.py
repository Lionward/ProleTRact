"""
FastAPI backend for ProleTRact React application
"""
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional, Dict, Any
import pysam
from pathlib import Path
import uvicorn
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
# pandas is only imported when we need it for the pathogenic catalog stuff

# how many workers to use for processing cohorts
# using all CPU cores since parsing VCF records is CPU intensive
COHORT_WORKERS = multiprocessing.cpu_count()

app = FastAPI(title="ProleTRact API")

# CORS stuff so the react frontend can talk to us
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", "http://localhost:8501"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# in-memory cache (should probably use redis or a db in production but this works for now)
vcf_cache = {}
# cache for cohort sample names
cohort_sample_cache = {}
# cache for cohort regions
cohort_regions_cache = {}

def parse_record_assembly(vcf_file: str, region: str) -> Optional[Dict[str, Any]]:
    """
    Parse an assembly VCF record (single haplotype per file).
    Used for population VCF files (h1/h2 files).
    """
    try:
        vcf = pysam.VariantFile(vcf_file)
        record_iter = vcf.fetch(region=region)
        rec = next(record_iter, None)
        vcf.close()
        
        if rec is None:
            return None
        
        # get motif ids for the ALT allele (this is a single haplotype file)
        ids_h = rec.samples[0].get("MI", None)
        if ids_h is None:
            ids_h_list = []
        elif isinstance(ids_h, (tuple, list)):
            ids_h_list = [str(x) for x in ids_h if x]
        else:
            ids_h_list = str(ids_h).split("_") if ids_h else []
        
        # get motif ids for the REF allele
        ids_ref = rec.info.get('MOTIF_IDs_REF', [])
        if isinstance(ids_ref, (tuple, list)):
            ids_ref = [str(x) for x in ids_ref if x]
        elif ids_ref:
            ids_ref = str(ids_ref).split("_")
        else:
            ids_ref = []
        
        # copy numbers for ref and alt alleles
        ref_CN = rec.info.get('CN_ref', 0)
        CN_H = rec.samples[0].get('CN', 0)
        if isinstance(CN_H, (tuple, list)):
            CN_H = CN_H[0] if len(CN_H) > 0 else 0
        
        # get motif names from INFO field
        motif_names = rec.info.get('MOTIFS', [])
        if isinstance(motif_names, tuple):
            motif_names = list(motif_names)
        elif not isinstance(motif_names, list):
            motif_names = [motif_names] if motif_names else []
        
        alt_allele = rec.alts[0] if rec.alts and rec.alts[0] != '.' else ''
        
        # get the span of the motifs
        spans = rec.samples[0].get('SP', "")
        if spans is None:
            spans = ""
        
        # get the genotype
        gt = rec.samples[0].get('GT', (0,))
        if isinstance(gt, (tuple, list)):
            gt_str = str(gt[0]) if len(gt) > 0 else "0"
        else:
            gt_str = str(gt)
        
        record = {
            'chr': rec.chrom,
            'pos': rec.pos,
            'stop': rec.stop,
            'motifs': motif_names,
            'motif_ids_h': ids_h_list,
            'motif_ids_ref': ids_ref,
            'ref_CN': ref_CN,
            'CN_H': CN_H,
            'spans': spans,
            'ref_allele': rec.ref,
            'alt_allele': alt_allele,
            'gt': gt_str,
            'id': rec.id,
        }
        return record
    except Exception as e:
        print(f"Error parsing assembly record: {e}")
        return None

def parse_record(vcf_file: str, region: str) -> Optional[Dict[str, Any]]:
    """
    Parse a single VCF record for a specified region.
    Similar to parsers.parse_record in the original ProleTRact.
    """
    try:
        vcf = pysam.VariantFile(vcf_file)
        record_iter = vcf.fetch(region=region)
        rec = next(record_iter, None)
        vcf.close()
        
        if rec is None:
            return None
        
        # parse motif IDs for h1 and h2
        mi = rec.samples[0].get('MI', None)
        if mi is None:
            ids_h1 = []
            ids_h2 = []
        elif isinstance(mi, tuple):
            ids_h1 = mi[0].split("_") if mi[0] else []
            ids_h2 = mi[1].split("_") if len(mi) > 1 and mi[1] else []
        else:
            ids_h1 = str(mi).split("_") if mi else []
            ids_h2 = ids_h1.copy()
        
        # figure out which alleles each haplotype has using the GT field
        ref_allele = rec.ref
        gt = rec.samples[0].get('GT', (0, 0))
        if not isinstance(gt, (tuple, list)):
            gt = (0, 0)
        
        # get all alleles - ref is 0, alt[0] is 1, alt[1] is 2, etc
        all_alleles = [ref_allele]
        if rec.alts:
            all_alleles.extend([alt for alt in rec.alts if alt != '.'])
        
        # figure out which allele each haplotype has from the GT
        gt_h1 = gt[0] if len(gt) > 0 else 0
        gt_h2 = gt[1] if len(gt) > 1 else gt_h1
        
        # get the actual sequences for each haplotype
        alt_allele1 = all_alleles[gt_h1] if gt_h1 < len(all_alleles) else ref_allele
        alt_allele2 = all_alleles[gt_h2] if gt_h2 < len(all_alleles) else ref_allele
        
        # Debug: Log sequence lengths for NA13509 sample
        if 'NA13509' in str(vcf_file):
            print(f"DEBUG NA13509: GT={gt}, all_alleles count={len(all_alleles)}")
            print(f"DEBUG NA13509: alt_allele1 length={len(alt_allele1) if alt_allele1 else 0}, alt_allele2 length={len(alt_allele2) if alt_allele2 else 0}")
            if len(all_alleles) > 1:
                print(f"DEBUG NA13509: ALT[0] length={len(all_alleles[1]) if len(all_alleles) > 1 else 0}")
            if len(all_alleles) > 2:
                print(f"DEBUG NA13509: ALT[1] length={len(all_alleles[2]) if len(all_alleles) > 2 else 0}")
        
        # removed debug logging, it was too slow
        
        # even if both haplotypes have the same allele sequence, we still show both
        # because they might have different motif IDs
        
        # copy numbers for h1 and h2
        CNs = rec.samples[0].get('CN', (0, 0))
        if isinstance(CNs, tuple):
            CN_H1 = str(CNs[0]) if len(CNs) > 0 else None
            CN_H2 = str(CNs[1]) if len(CNs) > 1 else None
        else:
            CN_H1 = str(CNs)
            CN_H2 = str(CNs)
        
        # parse the span info
        SP_field = rec.samples[0].get('SP', None)
        if SP_field is None:
            spans_h1 = ""
            spans_h2 = ""
        elif isinstance(SP_field, tuple):
            spans_h1 = SP_field[0] if len(SP_field) > 0 and SP_field[0] is not None else ""
            spans_h2 = SP_field[1] if len(SP_field) > 1 and SP_field[1] is not None else spans_h1
        else:
            spans_h1 = str(SP_field) if SP_field else ""
            spans_h2 = spans_h1
        
        # replace None with empty string
        ref_span = rec.info.get('REF_SPAN', None)
        spans = [
            str(ref_span) if ref_span is not None else "",
            spans_h1 if spans_h1 else "",
            spans_h2 if spans_h2 else ""
        ]
        
        # get motif names from INFO field
        motif_names = rec.info.get('MOTIFS', [])
        if isinstance(motif_names, tuple):
            motif_names = list(motif_names)
        elif not isinstance(motif_names, list):
            motif_names = [motif_names] if motif_names else []
        
        # get motif IDs for the ref allele
        motif_ids_ref = rec.info.get('MOTIF_IDs_REF', "")
        if motif_ids_ref:
            if isinstance(motif_ids_ref, (tuple, list)):
                motif_ids_ref = "_".join(str(x) for x in motif_ids_ref)
            motif_ids_ref = str(motif_ids_ref).split("_")
        else:
            motif_ids_ref = []
        
        # make a genotype string (we already got gt above)
        gt_str = '/'.join([str(i) for i in gt]) if isinstance(gt, (tuple, list)) else str(gt)
        
        # get supporting reads
        supporting_reads = rec.samples[0].get('DP', None)
        if supporting_reads is None:
            supporting_reads_h1 = 0
            supporting_reads_h2 = 0
        elif isinstance(supporting_reads, tuple):
            supporting_reads_h1 = supporting_reads[0] if len(supporting_reads) > 0 else 0
            supporting_reads_h2 = supporting_reads[1] if len(supporting_reads) > 1 else supporting_reads_h1
        else:
            supporting_reads_h1 = int(supporting_reads) if supporting_reads else 0
            supporting_reads_h2 = supporting_reads_h1
        
        # build the final record dict
        record = {
            'chr': rec.chrom,
            'pos': rec.pos,
            'stop': rec.stop,
            'motifs': motif_names,
            'motif_ids_h1': ids_h1,
            'motif_ids_h2': ids_h2,
            'motif_ids_ref': motif_ids_ref,
            'ref_CN': rec.info.get('CN_ref', None),
            'CN_H1': CN_H1,
            'CN_H2': CN_H2,
            'spans': spans,
            'ref_allele': ref_allele,
            'alt_allele1': alt_allele1,
            'alt_allele2': alt_allele2,
            'gt': gt_str,
            'supported_reads_h1': supporting_reads_h1,
            'supported_reads_h2': supporting_reads_h2,
            'id': rec.id,
        }
        
        # Debug: Verify sequences are not truncated for NA13509
        if 'NA13509' in str(vcf_file):
            print(f"DEBUG NA13509 record: alt_allele1 length={len(alt_allele1) if alt_allele1 else 0}, alt_allele2 length={len(alt_allele2) if alt_allele2 else 0}")
            print(f"DEBUG NA13509 record: alt_allele1 preview={alt_allele1[:50] if alt_allele1 else 'None'}...")
            print(f"DEBUG NA13509 record: alt_allele2 preview={alt_allele2[:50] if alt_allele2 else 'None'}...")
        
        return record
    except Exception as e:
        print(f"Error parsing record: {e}")
        return None

class VCFLoadRequest(BaseModel):
    vcf_path: str

class FilterRequest(BaseModel):
    vcf_path: str
    genotype_filter: Optional[List[str]] = None
    page: int = 0
    page_size: int = 50

class FilterAdvancedRequest(BaseModel):
    vcf_path: str
    motif_size_min: Optional[int] = None
    motif_size_max: Optional[int] = None
    cn_min: Optional[float] = None
    cn_max: Optional[float] = None
    chromosomes: Optional[List[str]] = None
    genotypes: Optional[List[str]] = None
    pathogenic_only: bool = False
    has_annotations: bool = False
    annotation_tags: Optional[List[str]] = None
    annotated_regions: Optional[List[str]] = None
    page: int = 0
    page_size: int = 50

class RegionInfo(BaseModel):
    id: str
    region: str
    genotype: str

class FilterResponse(BaseModel):
    records: List[RegionInfo]
    total_matching: int
    total_regions: int
    current_page: int
    total_pages: int
    available_genotypes: Optional[List[str]] = None

@app.get("/")
async def root():
    from proletract import __version__
    return {"message": "ProleTRact API", "version": __version__}


# --- File browser API (for modern file/folder selection) ---
import os

def _get_browse_roots() -> List[str]:
    """Return allowed root directories for browsing. Uses env var or defaults."""
    env_root = os.environ.get("PROLETRACT_BROWSE_ROOT")
    if env_root:
        return [str(Path(env_root).resolve())]
    roots = []
    if os.name != "nt":
        roots.append("/")
    home = os.path.expanduser("~")
    if home and home not in roots:
        roots.append(str(Path(home).resolve()))
    return roots or [str(Path.cwd())]


def _path_is_allowed(requested: Path) -> bool:
    """Check if path is within allowed browse roots (prevents path traversal)."""
    try:
        resolved = requested.resolve()
        for root in _get_browse_roots():
            root_path = Path(root).resolve()
            try:
                resolved.relative_to(root_path)
                return True
            except ValueError:
                continue
        return False
    except (OSError, RuntimeError):
        return False


@app.get("/api/files/roots")
async def get_file_browse_roots():
    """Return root directories available for file browsing."""
    return {"roots": [{"path": r, "name": Path(r).name or r} for r in _get_browse_roots()]}


@app.get("/api/files/browse")
async def browse_files(path: str = "", mode: str = "file"):
    """
    List directory contents for file browser.
    mode: 'file' = show VCF files + dirs (for selecting a VCF file)
          'folder' = show dirs only (for selecting a folder)
    """
    try:
        roots = _get_browse_roots()
        if not path:
            # Return roots as "directories" when no path specified
            dir_entries = []
            for r in roots:
                rp = Path(r)
                name = rp.name if rp.name else ("Home" if "home" in r.lower() else "Root")
                dir_entries.append({"name": name, "path": r})
            return {
                "path": "",
                "parent": None,
                "directories": dir_entries,
                "files": [],
            }
        p = Path(path)
        if not p.is_absolute():
            p = Path(roots[0]) / path if roots else Path.cwd() / path
        p = p.resolve()
        if not _path_is_allowed(p):
            raise HTTPException(status_code=403, detail="Path not allowed")
        if not p.exists():
            raise HTTPException(status_code=404, detail="Path not found")
        if not p.is_dir():
            raise HTTPException(status_code=400, detail="Not a directory")
        dirs = []
        files = []
        for item in sorted(p.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower())):
            try:
                if item.is_dir() and not item.name.startswith("."):
                    dirs.append({"name": item.name, "path": str(item)})
                elif item.is_file():
                    if mode == "file" and (item.name.endswith(".vcf") or item.name.endswith(".vcf.gz")):
                        files.append({"name": item.name, "path": str(item)})
                    elif mode == "folder":
                        files.append({"name": item.name, "path": str(item)})
            except (OSError, PermissionError):
                continue
        parent = str(p.parent) if str(p.parent) != str(p) else None
        if parent and not _path_is_allowed(Path(parent)):
            parent = None
        return {
            "path": str(p),
            "parent": parent,
            "directories": dirs,
            "files": files,
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/vcf/clear-cache")
async def clear_vcf_cache(vcf_path: Optional[str] = None):
    """Clear VCF cache for a specific file or all files"""
    try:
        if vcf_path:
            if vcf_path in vcf_cache:
                del vcf_cache[vcf_path]
                return {"success": True, "message": f"Cache cleared for {vcf_path}"}
            else:
                return {"success": False, "message": "VCF not found in cache"}
        else:
            vcf_cache.clear()
            return {"success": True, "message": "All VCF caches cleared"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/vcf/load")
async def load_vcf(request: VCFLoadRequest):
    """Load and parse VCF file - uses same approach as statistics to ensure consistency"""
    try:
        if not Path(request.vcf_path).exists():
            raise HTTPException(status_code=404, detail="VCF file not found")
        
        # clear cache for this vcf to make sure we get fresh data
        if request.vcf_path in vcf_cache:
            del vcf_cache[request.vcf_path]
        
        print(f"Loading VCF file: {request.vcf_path}")
        vcf = pysam.VariantFile(request.vcf_path)
        records = []
        region_genotypes = {}
        
        # use the same approach as the stats endpoint - fetch() without args
        # this way we get ALL records, same as what stats shows
        # stats can read 1.2M+ regions like this so we should be fine
        print("Reading all records from VCF file...")
        record_count = 0
        chromosomes_seen = set()
        for rec in vcf.fetch():
            record_count += 1
            region_str = f"{rec.chrom}:{rec.pos}-{rec.stop}"
            chromosomes_seen.add(rec.chrom)
            
            # motif size (max length of motifs in INFO)
            motifs = rec.info.get('MOTIFS', [])
            if isinstance(motifs, tuple):
                motifs = list(motifs)
            elif not isinstance(motifs, list):
                motifs = [motifs] if motifs else []
            motif_sizes = [len(str(m)) for m in motifs if m]
            motif_size = max(motif_sizes) if motif_sizes else 0
            
            # copy number (max of CN from samples, or CN_ref from INFO)
            cn_max_val = 0
            try:
                cn_ref = rec.info.get('CN_ref')
                if cn_ref is not None:
                    cn_max_val = max(cn_max_val, float(cn_ref))
                cn = rec.samples[0].get('CN')
                if cn is not None:
                    if isinstance(cn, (tuple, list)):
                        for c in cn:
                            if c is not None:
                                cn_max_val = max(cn_max_val, float(c))
                    else:
                        cn_max_val = max(cn_max_val, float(cn))
            except (TypeError, ValueError, KeyError, IndexError):
                pass
            
            records.append({
                'id': rec.id,
                'region': region_str,
                'chrom': rec.chrom,
                'pos': rec.pos,
                'stop': rec.stop,
                'motif_size': motif_size,
                'cn_max': cn_max_val
            })
            
            # extract the genotype
            try:
                gt = rec.samples[0]['GT']
                if gt is not None:
                    gt_str = '/'.join([str(i) for i in gt]) if isinstance(gt, (tuple, list)) else str(gt)
                else:
                    gt_str = './.'
                region_genotypes[region_str] = gt_str
            except (KeyError, IndexError, AttributeError):
                region_genotypes[region_str] = './.'
            
            # print progress every 100k records
            if record_count % 100000 == 0:
                print(f"  Loaded {record_count:,} regions...")
        
        vcf.close()
        
        print(f"Total regions loaded: {len(records):,}")
        
        # cache the results
        vcf_cache[request.vcf_path] = {
            'records': records,
            'region_genotypes': region_genotypes
        }
        
        available_genotypes = sorted(set(region_genotypes.values()))
        # Natural sort: chr1, chr2, ..., chr22, chrX, chrY, chrM
        def chrom_sort_key(c):
            s = str(c).replace('chr', '').replace('Chr', '')
            if s == 'X': return (23, 'X')
            if s == 'Y': return (24, 'Y')
            if s in ('M', 'MT'): return (25, s)
            try: return (int(s), s)
            except ValueError: return (999, s)
        available_chromosomes = sorted(chromosomes_seen, key=chrom_sort_key)
        
        return {
            "success": True,
            "total_regions": len(records),
            "available_genotypes": available_genotypes,
            "available_chromosomes": available_chromosomes,
            "message": f"Loaded {len(records):,} regions"
        }
    except Exception as e:
        print(f"Error loading VCF: {e}")
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/vcf/filter", response_model=FilterResponse)
async def filter_regions(request: FilterRequest):
    """Filter regions with server-side pagination"""
    try:
        if request.vcf_path not in vcf_cache:
            raise HTTPException(status_code=404, detail="VCF not loaded. Load VCF first.")
        
        cache = vcf_cache[request.vcf_path]
        records = cache['records']
        region_genotypes = cache['region_genotypes']
        
        # apply genotype filter if specified
        if request.genotype_filter:
            filtered = [
                r for r in records
                if region_genotypes.get(r['region'], './.') in request.genotype_filter
            ]
        else:
            filtered = records
        
        total_matching = len(filtered)
        
        # pagination
        start_idx = request.page * request.page_size
        end_idx = start_idx + request.page_size
        page_records = filtered[start_idx:end_idx]
        
        # format the response
        result_records = [
            RegionInfo(
                id=r['id'],
                region=r['region'],
                genotype=region_genotypes.get(r['region'], './.')
            )
            for r in page_records
        ]
        
        total_pages = (total_matching // request.page_size) + (1 if total_matching % request.page_size > 0 else 0)
        available_genotypes = sorted(set(region_genotypes.values()))
        
        return FilterResponse(
            records=result_records,
            total_matching=total_matching,
            total_regions=len(records),
            current_page=request.page,
            total_pages=total_pages,
            available_genotypes=available_genotypes
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

def _is_pathogenic(region_str: str, cn_max: float, catalog) -> bool:
    """Check if region is pathogenic (overlaps catalog and CN >= threshold)."""
    if catalog is None:
        return False
    try:
        parts = region_str.split(':')
        if len(parts) < 2:
            return False
        chrom = parts[0]
        pos_end = parts[1].split('-')
        if len(pos_end) < 2:
            return False
        start = int(pos_end[0])
        end = int(pos_end[1])
        chr_norm = chrom if chrom.startswith('chr') else f'chr{chrom}'
        if isinstance(catalog, list):
            for entry in catalog:
                if entry.get('chrom') not in (chrom, chr_norm):
                    continue
                if (entry['start'] <= end) and (entry['end'] >= start):
                    thresh = entry.get('pathogenic_min')
                    if thresh is not None and cn_max >= float(thresh):
                        return True
            return False
        import pandas as pd
        if isinstance(catalog, pd.DataFrame) and not catalog.empty:
            chr_matches = catalog[(catalog['chrom'] == chrom) | (catalog['chrom'] == chr_norm)]
            for _, row in chr_matches.iterrows():
                if (row['start'] <= end) and (row['end'] >= start):
                    thresh = row.get('pathogenic_min')
                    if pd.notna(thresh) and cn_max >= float(thresh):
                        return True
        return False
    except Exception:
        return False

@app.post("/api/vcf/filter-advanced", response_model=FilterResponse)
async def filter_regions_advanced(request: FilterAdvancedRequest):
    """Filter regions with advanced criteria (motif size, CN, chromosomes, genotypes, pathogenic)."""
    try:
        if request.vcf_path not in vcf_cache:
            raise HTTPException(status_code=404, detail="VCF not loaded. Load VCF first.")
        
        cache = vcf_cache[request.vcf_path]
        records = cache['records']
        region_genotypes = cache['region_genotypes']
        
        # Build filtered list
        filtered = []
        pathogenic_catalog = load_pathogenic_catalog() if request.pathogenic_only else None
        
        for r in records:
            region_str = r['region']
            gt = region_genotypes.get(region_str, './.')
            
            # genotype filter
            if request.genotypes and len(request.genotypes) > 0:
                if gt not in request.genotypes:
                    continue
            
            # chromosome filter
            if request.chromosomes and len(request.chromosomes) > 0:
                chrom = r.get('chrom', '')
                if chrom not in request.chromosomes:
                    continue
            
            # motif size filter
            motif_size = r.get('motif_size', 0)
            if request.motif_size_min is not None and motif_size < request.motif_size_min:
                continue
            if request.motif_size_max is not None and motif_size > request.motif_size_max:
                continue
            
            # copy number filter
            cn_max = r.get('cn_max', 0)
            if request.cn_min is not None and cn_max < request.cn_min:
                continue
            if request.cn_max is not None and cn_max > request.cn_max:
                continue
            
            # pathogenic filter
            if request.pathogenic_only:
                if not _is_pathogenic(region_str, cn_max, pathogenic_catalog):
                    continue
            
            # annotated_regions filter (frontend sends client-filtered list)
            if request.annotated_regions and len(request.annotated_regions) > 0:
                if region_str not in request.annotated_regions:
                    continue
            
            filtered.append(r)
        
        total_matching = len(filtered)
        start_idx = request.page * request.page_size
        end_idx = start_idx + request.page_size
        page_records = filtered[start_idx:end_idx]
        
        result_records = [
            RegionInfo(
                id=r['id'],
                region=r['region'],
                genotype=region_genotypes.get(r['region'], './.')
            )
            for r in page_records
        ]
        
        total_pages = (total_matching // request.page_size) + (1 if total_matching % request.page_size > 0 else 0)
        available_genotypes = sorted(set(region_genotypes.values()))
        
        return FilterResponse(
            records=result_records,
            total_matching=total_matching,
            total_regions=len(records),
            current_page=request.page,
            total_pages=total_pages,
            available_genotypes=available_genotypes
        )
    except HTTPException:
        raise
    except Exception as e:
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/vcf/regions")
async def get_all_regions(vcf_path: str):
    """Get all available regions for autocomplete"""
    try:
        if not Path(vcf_path).exists():
            raise HTTPException(status_code=404, detail="VCF file not found")
        
        if vcf_path not in vcf_cache:
            raise HTTPException(status_code=404, detail="VCF not loaded. Load VCF first.")
        
        cache = vcf_cache[vcf_path]
        regions = [r['region'] for r in cache['records']]
        
        return {
            "success": True,
            "regions": regions
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/vcf/region-page")
async def get_region_page(vcf_path: str, region: str, genotype_filter: Optional[str] = None, page_size: int = 50):
    """Find which page a specific region is on"""
    try:
        if not Path(vcf_path).exists():
            raise HTTPException(status_code=404, detail="VCF file not found")
        
        if vcf_path not in vcf_cache:
            raise HTTPException(status_code=404, detail="VCF not loaded. Load VCF first.")
        
        cache = vcf_cache[vcf_path]
        records = cache['records']
        region_genotypes = cache['region_genotypes']
        
        # Apply genotype filter
        genotype_list = None
        if genotype_filter:
            genotype_list = genotype_filter.split(',')
            filtered = [
                r for r in records
                if region_genotypes.get(r['region'], './.') in genotype_list
            ]
        else:
            filtered = records
        
        # Find the index of the region
        region_index = None
        for idx, r in enumerate(filtered):
            if r['region'] == region:
                region_index = idx
                break
        
        if region_index is None:
            raise HTTPException(status_code=404, detail="Region not found in filtered results")
        
        # Calculate page number (0-indexed)
        page_number = region_index // page_size
        
        return {
            "success": True,
            "page": page_number,
            "index": region_index,
            "total_matching": len(filtered)
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/vcf/region-by-index")
async def get_region_by_index(vcf_path: str, region_index: int, genotype_filter: Optional[str] = None, page_size: int = 50):
    """Get a region by its index in the filtered list"""
    try:
        if not Path(vcf_path).exists():
            raise HTTPException(status_code=404, detail="VCF file not found")
        
        if vcf_path not in vcf_cache:
            raise HTTPException(status_code=404, detail="VCF not loaded. Load VCF first.")
        
        cache = vcf_cache[vcf_path]
        records = cache['records']
        region_genotypes = cache['region_genotypes']
        
        # Apply genotype filter
        if genotype_filter:
            genotype_list = genotype_filter.split(',')
            filtered = [
                r for r in records
                if region_genotypes.get(r['region'], './.') in genotype_list
            ]
        else:
            filtered = records
        
        # Validate index
        if region_index < 0 or region_index >= len(filtered):
            raise HTTPException(status_code=404, detail=f"Region index {region_index} out of range (0-{len(filtered)-1})")
        
        # Get the region at this index
        target_region = filtered[region_index]
        
        # Calculate page number (0-indexed)
        page_number = region_index // page_size
        
        return {
            "success": True,
            "region": target_region['region'],
            "page": page_number,
            "index": region_index,
            "total_matching": len(filtered)
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/vcf/region/{region_str}")
async def get_region_data(region_str: str, vcf_path: str):
    """Get detailed data for a specific region"""
    try:
        if not Path(vcf_path).exists():
            raise HTTPException(status_code=404, detail="VCF file not found")
        
        record = parse_record(vcf_path, region_str)
        if record is None:
            raise HTTPException(status_code=404, detail="Region not found")
        
        return {
            "success": True,
            "record": record
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/vcf/statistics")
async def get_vcf_statistics(vcf_path: str):
    """Get comprehensive statistics about the VCF file"""
    try:
        if not Path(vcf_path).exists():
            raise HTTPException(status_code=404, detail="VCF file not found")
        
        vcf = pysam.VariantFile(vcf_path)
        
        # Collect statistics
        total_regions = 0
        motif_size_counts = {}  # max_motif_size -> count
        motif_lengths = []  # all motif lengths for distribution
        regions_by_chromosome = {}
        genotype_counts = {}  # genotype -> count
        
        # Get all chromosomes from the VCF header to ensure we read all records
        # This is important for large indexed VCF files where fetch() without args might miss records
        chromosomes = list(vcf.header.contigs.keys())
        
        # If we have chromosomes, iterate through each one explicitly
        # This ensures we read ALL records, even in very large indexed VCF files
        if chromosomes:
            for chrom in chromosomes:
                try:
                    for rec in vcf.fetch(chrom):
                        total_regions += 1
                        
                        # Get motifs from INFO field
                        motifs = rec.info.get('MOTIFS', [])
                        if isinstance(motifs, tuple):
                            motifs = list(motifs)
                        elif not isinstance(motifs, list):
                            motifs = [motifs] if motifs else []
                        
                        # Calculate max motif size for this region
                        max_motif_size = 0
                        if motifs:
                            motif_lengths_in_region = [len(str(m)) for m in motifs if m]
                            if motif_lengths_in_region:
                                max_motif_size = max(motif_lengths_in_region)
                                motif_lengths.extend(motif_lengths_in_region)
                        
                        # Categorize motif size
                        if max_motif_size == 0:
                            category = "Unknown"
                        elif max_motif_size <= 10:
                            category = str(max_motif_size)
                        else:
                            category = ">10"
                        motif_size_counts[category] = motif_size_counts.get(category, 0) + 1
                        
                        # Count by chromosome
                        chrom_name = rec.chrom
                        regions_by_chromosome[chrom_name] = regions_by_chromosome.get(chrom_name, 0) + 1
                        
                        # Extract genotype information
                        if len(rec.samples) > 0 and 'GT' in rec.samples[0]:
                            gt = rec.samples[0]['GT']
                            if gt is not None:
                                if isinstance(gt, tuple):
                                    gt_str = '/'.join([str(g) if g is not None else '.' for g in gt])
                                elif isinstance(gt, list):
                                    gt_str = '/'.join([str(g) if g is not None else '.' for g in gt])
                                else:
                                    gt_str = str(gt)
                                genotype_counts[gt_str] = genotype_counts.get(gt_str, 0) + 1
                except (ValueError, KeyError):
                    # Skip chromosomes that don't exist or can't be fetched
                    continue
        else:
            # Fallback: if no chromosomes in header, use fetch() without arguments
            for rec in vcf.fetch():
                total_regions += 1
                
                # Get motifs from INFO field
                motifs = rec.info.get('MOTIFS', [])
                if isinstance(motifs, tuple):
                    motifs = list(motifs)
                elif not isinstance(motifs, list):
                    motifs = [motifs] if motifs else []
                
                # Calculate max motif size for this region
                max_motif_size = 0
                if motifs:
                    motif_lengths_in_region = [len(str(m)) for m in motifs if m]
                    if motif_lengths_in_region:
                        max_motif_size = max(motif_lengths_in_region)
                        motif_lengths.extend(motif_lengths_in_region)
                
                # Categorize motif size
                if max_motif_size == 0:
                    category = "Unknown"
                elif max_motif_size <= 10:
                    category = str(max_motif_size)
                else:
                    category = ">10"
                motif_size_counts[category] = motif_size_counts.get(category, 0) + 1
                
                # Count by chromosome
                chrom = rec.chrom
                regions_by_chromosome[chrom] = regions_by_chromosome.get(chrom, 0) + 1
                
                # Extract genotype information
                if len(rec.samples) > 0 and 'GT' in rec.samples[0]:
                    gt = rec.samples[0]['GT']
                    if gt is not None:
                        if isinstance(gt, tuple):
                            gt_str = '/'.join([str(g) if g is not None else '.' for g in gt])
                        elif isinstance(gt, list):
                            gt_str = '/'.join([str(g) if g is not None else '.' for g in gt])
                        else:
                            gt_str = str(gt)
                        genotype_counts[gt_str] = genotype_counts.get(gt_str, 0) + 1
        
        vcf.close()
        
        # Calculate average and max motif sizes
        avg_motif_size = sum(motif_lengths) / len(motif_lengths) if motif_lengths else 0
        max_overall_motif_size = max(motif_lengths) if motif_lengths else 0
        
        # Bin motif lengths for histogram (to avoid sending huge arrays to frontend)
        motif_length_histogram = {}
        if motif_lengths:
            min_len = min(motif_lengths)
            max_len = max(motif_lengths)
            
            if min_len == max_len:
                # All values are the same
                motif_length_histogram = {str(int(min_len)): len(motif_lengths)}
            else:
                # Determine number of bins (between 15 and 30 for good histogram visualization)
                num_bins = min(30, max(15, int(len(motif_lengths) ** 0.5)))
                
                # Round min/max to nice numbers for cleaner bins
                # Use floor for min, ceiling for max to ensure all data is included
                min_rounded = int(min_len)
                max_rounded = int(max_len) + 1
                range_size = max_rounded - min_rounded
                
                # Calculate bin width (ensure it's at least 1)
                bin_width = max(1.0, range_size / num_bins)
                # Round bin width to a nice number
                if bin_width < 5:
                    bin_width = round(bin_width)
                else:
                    bin_width = round(bin_width / 5) * 5
                
                # Create bins and count values
                bins_dict = {}
                for length in motif_lengths:
                    # Calculate which bin this value belongs to
                    bin_idx = int((length - min_rounded) / bin_width)
                    bin_idx = min(bin_idx, num_bins - 1)  # Ensure it's within range
                    
                    # Calculate actual bin boundaries
                    bin_start = min_rounded + (bin_idx * bin_width)
                    bin_end = bin_start + bin_width
                    
                    # For the last bin, extend to include max value
                    if bin_idx == num_bins - 1:
                        bin_end = max_rounded
                    
                    # Create clean label
                    bin_start_int = int(bin_start)
                    bin_end_int = int(bin_end)
                    
                    if bin_end_int - bin_start_int <= 1:
                        bin_label = str(bin_start_int)
                    else:
                        bin_label = f"{bin_start_int}-{bin_end_int-1}"
                    
                    bins_dict[bin_label] = bins_dict.get(bin_label, 0) + 1
                
                # Only include bins with data (for cleaner visualization)
                motif_length_histogram = {k: v for k, v in bins_dict.items() if v > 0}
        
        return {
            "success": True,
            "total_regions": total_regions,
            "num_chromosomes": len(regions_by_chromosome),
            "avg_motif_size": round(avg_motif_size, 1),
            "max_motif_size": max_overall_motif_size,
            "motif_size_counts": motif_size_counts,
            "motif_length_histogram": motif_length_histogram,  # Changed from motif_lengths array to binned histogram
            "regions_by_chromosome": regions_by_chromosome,
            "genotype_counts": genotype_counts
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

class PopulationLoadRequest(BaseModel):
    folder_path: str

def process_vcf_file_for_loading(args):
    """Helper function to process a single VCF file for loading - used for parallel processing"""
    file_path_str = args
    try:
        file_path = Path(file_path_str)
        # Open VCF to extract sample name from header
        vcf = pysam.VariantFile(str(file_path))
        samples = list(vcf.header.samples)
        vcf.close()
        
        # Use sample name from VCF header, or filename if no samples
        if samples:
            sample_name = samples[0]  # Use first sample name
        else:
            # Fallback to filename without extension
            sample_name = file_path.stem.replace('.vcf', '')
        
        # Check if this is a haplotype-specific file
        is_haplotype, base_name, hap_suffix = is_haplotype_specific_name(sample_name)
        
        return {
            'filename': file_path.name,
            'sample_name': sample_name,
            'base_sample_name': base_name if is_haplotype else sample_name,
            'is_haplotype': is_haplotype,
            'haplotype_suffix': hap_suffix if is_haplotype else '',
            'path': str(file_path)
        }
    except Exception as e:
        # Skip files that can't be opened
        print(f"Warning: Could not read VCF header from {file_path_str}: {e}")
        return None

@app.post("/api/population/load")
async def load_population_vcf_files(request: PopulationLoadRequest):
    """Load population/cohort VCF files from a folder with parallel processing"""
    try:
        folder_path = Path(request.folder_path)
        if not folder_path.exists() or not folder_path.is_dir():
            raise HTTPException(status_code=404, detail="Folder not found or is not a directory")
        
        # Find all VCF files (.vcf.gz or .vcf)
        all_vcf_files = [f for f in folder_path.iterdir() 
                        if f.is_file() and (f.name.endswith('.vcf.gz') or f.name.endswith('.vcf'))]
        
        if not all_vcf_files:
            raise HTTPException(status_code=404, detail="No VCF files found in the folder")
        
        vcf_files = []
        sample_info = []
        
        # Process files in parallel using ProcessPoolExecutor for true parallelism
        file_paths = [str(f) for f in all_vcf_files]
        
        with ProcessPoolExecutor(max_workers=COHORT_WORKERS) as executor:
            # Submit all tasks
            future_to_file = {executor.submit(process_vcf_file_for_loading, file_path): file_path 
                            for file_path in file_paths}
            
            # Collect results as they complete
            for future in as_completed(future_to_file):
                try:
                    result = future.result()
                    if result is not None:
                        vcf_files.append(result['filename'])
                        sample_info.append(result)
                except Exception as e:
                    file_path = future_to_file[future]
                    print(f"Error processing {file_path}: {e}")
                    continue
        
        # Cache sample info for fast access later
        cohort_sample_cache[str(folder_path)] = sample_info
        
        return {
            "success": True,
            "file_count": len(vcf_files),
            "files": vcf_files,
            "sample_info": sample_info,
            "folder_path": str(folder_path)
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

def extract_regions_from_vcf_file(args):
    """Helper function to extract regions from a single VCF file - used for parallel processing"""
    vcf_file_path = args
    try:
        regions = set()
        vcf = pysam.VariantFile(str(vcf_file_path))
        for rec in vcf.fetch():
            region_str = f"{rec.chrom}:{rec.pos}-{rec.stop}"
            regions.add(region_str)
        vcf.close()
        return list(regions)
    except Exception as e:
        print(f"Warning: Could not read regions from {vcf_file_path}: {e}")
        return []

@app.get("/api/population/regions")
async def get_population_regions(folder_path: str):
    """Get all available regions from all VCF files in a cohort folder for autocomplete"""
    try:
        folder = Path(folder_path)
        if not folder.exists():
            raise HTTPException(status_code=404, detail="Population folder not found")
        
        # Check cache first
        cache_key = str(folder_path)
        if cache_key in cohort_regions_cache:
            cached_regions = cohort_regions_cache[cache_key]
            return {
                "success": True,
                "regions": cached_regions,
                "count": len(cached_regions),
                "cached": True
            }
        
        # Find all VCF files in the folder (.vcf.gz or .vcf)
        vcf_files = list(folder.glob("*.vcf.gz")) + list(folder.glob("*.vcf"))
        
        if not vcf_files:
            raise HTTPException(status_code=404, detail="No VCF files found in the folder")
        
        # Process files in parallel - MUCH faster than sequential
        all_regions = set()
        file_paths = [str(f) for f in vcf_files]
        
        # Use ProcessPoolExecutor for parallel processing
        with ProcessPoolExecutor(max_workers=COHORT_WORKERS) as executor:
            future_to_file = {executor.submit(extract_regions_from_vcf_file, file_path): file_path 
                            for file_path in file_paths}
            
            # Collect results as they complete
            for future in as_completed(future_to_file):
                try:
                    regions = future.result()
                    all_regions.update(regions)
                except Exception as e:
                    file_path = future_to_file[future]
                    print(f"Error extracting regions from {file_path}: {e}")
                    continue
        
        # Sort regions properly by chromosome and position
        def sort_region(region_str):
            """Parse region string and return sortable tuple"""
            match = re.match(r'^([^:]+):(\d+)-(\d+)$', region_str)
            if not match:
                return (999, '', 0, 0)
            chr_name = match.group(1)
            start = int(match.group(2))
            end = int(match.group(3))
            # Extract numeric part of chromosome for proper sorting
            chr_num_str = chr_name.replace('chr', '').replace('Chr', '').replace('CHR', '')
            try:
                chr_num = int(chr_num_str)
            except ValueError:
                # Handle X, Y, M, etc.
                chr_map = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}
                chr_num = chr_map.get(chr_num_str.upper(), 999)
            return (chr_num, chr_name, start, end)
        
        sorted_regions = sorted(list(all_regions), key=sort_region)
        
        # Cache the results
        cohort_regions_cache[cache_key] = sorted_regions
        
        return {
            "success": True,
            "regions": sorted_regions,
            "count": len(sorted_regions),
            "cached": False
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

def get_sample_id_only(args):
    """Helper function to get just sample ID and basic metadata - used for fast initial loading"""
    vcf_file_path, region_str = args
    # REMOVED: Debug logging was causing significant performance issues
    try:
        vcf_file = Path(vcf_file_path)
        # Extract sample name from VCF header
        vcf = pysam.VariantFile(str(vcf_file))
        samples = list(vcf.header.samples)
        
        # Check if region exists (quick check without full parsing)
        record_iter = vcf.fetch(region=region_str)
        rec = next(record_iter, None)
        vcf.close()
        
        if rec is None:
            return None
        
        # Use sample name from VCF header, or filename if no samples
        if samples:
            sample_name = samples[0]  # Use first sample name
        else:
            # Fallback to filename without extension
            sample_name = vcf_file.stem.replace('.vcf', '')
        
        # Return minimal metadata
        result = {
            'sample_name': sample_name,
            'file_path': str(vcf_file),
            'has_data': True
        }
        return result
    except Exception as e:
        # Skip files that can't be parsed or don't have the region
        return None

def is_haplotype_specific_name(name: str):
    """
    Check if a sample name indicates a haplotype-specific file.
    Returns: (is_haplotype_specific, base_name, haplotype_suffix)
    Examples:
    - 'sample_h1' -> (True, 'sample', '_h1')
    - 'sample_hap1' -> (True, 'sample', '_hap1')
    - 'sample_haplotype_1' -> (True, 'sample', '_haplotype_1')
    - 'sample' -> (False, 'sample', '')
    """
    import re
    # Patterns for haplotype-specific naming
    patterns = [
        (r'^(.+?)_h([12])$', r'_h\2'),  # sample_h1, sample_h2
        (r'^(.+?)_hap([12])$', r'_hap\2'),  # sample_hap1, sample_hap2
        (r'^(.+?)_haplotype[_-]?([12])$', r'_haplotype_\2'),  # sample_haplotype_1, sample_haplotype-1
    ]
    
    for pattern, suffix_pattern in patterns:
        match = re.match(pattern, name, re.IGNORECASE)
        if match:
            base_name = match.group(1)
            hap_num = match.group(2)
            suffix = suffix_pattern.replace(r'\2', hap_num)
            return (True, base_name, suffix)
    
    return (False, name, '')

def process_single_vcf_file(args):
    """Helper function to process a single VCF file - used for parallel processing"""
    vcf_file_path, region_str, cohort_mode = args
    try:
        vcf_file = Path(vcf_file_path)
        # Extract sample name from VCF header
        vcf = pysam.VariantFile(str(vcf_file))
        samples = list(vcf.header.samples)
        vcf.close()
        
        # Use sample name from VCF header, or filename if no samples
        if samples:
            sample_name = samples[0]  # Use first sample name
        else:
            # Fallback to filename without extension
            sample_name = vcf_file.stem.replace('.vcf', '')
        
        # Parse region to check chromosome
        region_parts = region_str.split(':')
        chrom = region_parts[0] if len(region_parts) > 0 else ''
        is_sex_chrom = chrom.upper() in ['X', 'Y', 'CHRX', 'CHRY']
        
        # Use the mode parameter to determine parsing strategy
        # cohort_mode can be 'cohort-read' (diploid) or 'cohort-assembly' (single haplotype)
        record = None
        if cohort_mode == 'cohort-read':
            # Read-based mode: use diploid parsing (parse_record)
            record = parse_record(str(vcf_file), region_str)
        elif cohort_mode == 'cohort-assembly':
            # Assembly-based mode: use assembly parsing (parse_record_assembly)
            record = parse_record_assembly(str(vcf_file), region_str)
        else:
            # Fallback: auto-detect based on GT field (for backward compatibility)
            try:
                vcf_check = pysam.VariantFile(str(vcf_file))
                record_iter = vcf_check.fetch(region=region_str)
                rec_check = next(record_iter, None)
                vcf_check.close()
                
                if rec_check:
                    gt_check = rec_check.samples[0].get('GT', None)
                    is_diploid_gt = isinstance(gt_check, (tuple, list)) and len(gt_check) == 2
                    
                    if is_sex_chrom:
                        record = parse_record_assembly(str(vcf_file), region_str)
                    elif is_diploid_gt:
                        record = parse_record(str(vcf_file), region_str)
                    else:
                        record = parse_record_assembly(str(vcf_file), region_str)
                else:
                    if is_sex_chrom:
                        record = parse_record_assembly(str(vcf_file), region_str)
                    else:
                        record = parse_record(str(vcf_file), region_str)
                        if not record:
                            record = parse_record_assembly(str(vcf_file), region_str)
            except Exception as e:
                # Fallback: try diploid first, then assembly
                if is_sex_chrom:
                    record = parse_record_assembly(str(vcf_file), region_str)
                else:
                    record = parse_record(str(vcf_file), region_str)
                    if not record:
                        record = parse_record_assembly(str(vcf_file), region_str)
        
        if record:
            # REMOVED: Debug logging was causing significant performance issues
            # Return with sample name
            return (sample_name, record)
        return None
    except Exception as e:
        # Skip files that can't be parsed or don't have the region
        print(f"Warning: Could not parse {vcf_file_path} for region {region_str}: {e}")
        return None

def get_sample_names_from_one_file(folder_path: str):
    """Get all sample names from one VCF file in the folder (fast - just reads headers)"""
    try:
        folder = Path(folder_path)
        vcf_files = list(folder.glob("*.vcf.gz")) + list(folder.glob("*.vcf"))
        
        if not vcf_files:
            return []
        
        # Just read the first file to get sample names
        first_file = vcf_files[0]
        vcf = pysam.VariantFile(str(first_file))
        samples = list(vcf.header.samples)
        vcf.close()
        
        # If no samples in header, use filename
        if not samples:
            samples = [first_file.stem.replace('.vcf', '')]
        
        # Get all file paths for mapping
        sample_to_file = {}
        for vcf_file in vcf_files:
            try:
                vcf = pysam.VariantFile(str(vcf_file))
                file_samples = list(vcf.header.samples)
                vcf.close()
                
                sample_name = file_samples[0] if file_samples else vcf_file.stem.replace('.vcf', '')
                sample_to_file[sample_name] = str(vcf_file)
            except:
                continue
        
        # Return sample names with their file paths
        return [{'sample_name': name, 'file_path': sample_to_file.get(name, '')} 
                for name in samples if name in sample_to_file]
    except Exception as e:
        print(f"Error getting sample names: {e}")
        return []

def get_sample_name_from_file(vcf_file_path: str):
    """Helper to get sample name from a single VCF file header (for parallel processing)"""
    try:
        vcf_file = Path(vcf_file_path)
        vcf = pysam.VariantFile(str(vcf_file))
        samples = list(vcf.header.samples)
        vcf.close()
        
        sample_name = samples[0] if samples else vcf_file.stem.replace('.vcf', '')
        return {
            'sample_name': sample_name,
            'file_path': str(vcf_file)
        }
    except Exception as e:
        print(f"Warning: Could not read header from {vcf_file_path}: {e}")
        return None

@app.get("/api/population/region/{region_str}/ids")
async def get_population_region_ids(region_str: str, folder_path: str, mode: str = 'cohort-read'):
    """Get sample IDs by loading one sample first, then return all sample names from folder"""
    import time
    start_time = time.time()
    
    try:
        folder = Path(folder_path)
        if not folder.exists():
            raise HTTPException(status_code=404, detail="Population folder not found")
        
        # Try to use cached sample info first (from /api/population/load)
        sample_ids = []
        if str(folder_path) in cohort_sample_cache:
            # Use cached sample info (much faster!)
            cached_info = cohort_sample_cache[str(folder_path)]
            sample_ids = [{'sample_name': info['sample_name'], 'file_path': info['path']} 
                         for info in cached_info]
            print(f"Using cached sample info: {len(sample_ids)} samples")
        else:
            # Fallback: find all VCF files and read headers in parallel
            vcf_files = list(folder.glob("*.vcf.gz")) + list(folder.glob("*.vcf"))
            
            if not vcf_files:
                raise HTTPException(status_code=404, detail="No VCF files found in the folder")
            
            file_paths = [str(f) for f in vcf_files]
            
            # Use ProcessPoolExecutor to read headers in parallel
            with ProcessPoolExecutor(max_workers=COHORT_WORKERS) as executor:
                future_to_file = {executor.submit(get_sample_name_from_file, file_path): file_path 
                                for file_path in file_paths}
                
                for future in as_completed(future_to_file):
                    try:
                        result = future.result()
                        if result is not None:
                            sample_ids.append(result)
                    except Exception as e:
                        file_path = future_to_file[future]
                        print(f"Error getting sample name from {file_path}: {e}")
                        continue
        
        if not sample_ids:
            raise HTTPException(status_code=404, detail="No samples found")
        
        # Load one sample first to get region data
        first_file_path = sample_ids[0]['file_path']
        first_record = None
        first_sample_name = sample_ids[0]['sample_name']
        
        first_sample_start = time.time()
        try:
            # Use the mode parameter to determine parsing
            if mode == 'cohort-read':
                record = parse_record(first_file_path, region_str)
            elif mode == 'cohort-assembly':
                record = parse_record_assembly(first_file_path, region_str)
            else:
                # Fallback: try both
                record = parse_record_assembly(first_file_path, region_str)
                if not record:
                    record = parse_record(first_file_path, region_str)
            
            if record:
                first_record = record
        except Exception as e:
            print(f"Warning: Could not load first sample from {first_file_path}: {e}")
        first_sample_time = time.time() - first_sample_start
        
        total_time = time.time() - start_time
        print(f"Timing - First sample: {first_sample_time:.2f}s, Total: {total_time:.2f}s, Samples: {len(sample_ids)}")
        
        return {
            "success": True,
            "sample_ids": sample_ids,
            "total": len(sample_ids),
            "first_sample": first_sample_name,
            "first_record": first_record  # Include first sample's data if available
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/population/region/{region_str}/samples")
async def get_population_region_samples(region_str: str, folder_path: str, sample_names: str = "", mode: str = 'cohort-read'):
    """Get full records for specific sample names (lazy loading) using multiprocessing"""
    import time
    start_time = time.time()
    
    try:
        folder = Path(folder_path)
        if not folder.exists():
            raise HTTPException(status_code=404, detail="Population folder not found")
        
        # Parse sample names from query parameter (comma-separated)
        requested_samples = [s.strip() for s in sample_names.split(',') if s.strip()] if sample_names else []
        
        if not requested_samples:
            raise HTTPException(status_code=400, detail="No sample names provided")
        
        # Use cached sample info if available (much faster!)
        sample_to_file = {}
        if str(folder_path) in cohort_sample_cache:
            cached_info = cohort_sample_cache[str(folder_path)]
            for info in cached_info:
                if info['sample_name'] in requested_samples:
                    sample_to_file[info['sample_name']] = info['path']
        else:
            # Fallback: scan files (slow)
            vcf_files = list(folder.glob("*.vcf.gz")) + list(folder.glob("*.vcf"))
            for vcf_file in vcf_files:
                try:
                    vcf = pysam.VariantFile(str(vcf_file))
                    samples = list(vcf.header.samples)
                    vcf.close()
                    
                    sample_name = samples[0] if samples else vcf_file.stem.replace('.vcf', '')
                    if sample_name in requested_samples:
                        sample_to_file[sample_name] = str(vcf_file)
                except Exception:
                    continue
        
        if not sample_to_file:
            return {
                "success": True,
                "records": {}
            }
        
        population_records = {}
        
        # Process only requested samples using ProcessPoolExecutor for true parallelism
        file_args = [(sample_to_file[sample], region_str, mode) for sample in requested_samples if sample in sample_to_file]
        
        # Use ProcessPoolExecutor for CPU-bound parsing tasks
        with ProcessPoolExecutor(max_workers=COHORT_WORKERS) as executor:
            future_to_file = {executor.submit(process_single_vcf_file, args): args[0] for args in file_args}
            
            for future in as_completed(future_to_file):
                try:
                    result = future.result()
                    if result is not None:
                        sample_name, record = result
                        population_records[sample_name] = record
                except Exception as e:
                    file_path = future_to_file[future]
                    print(f"Error processing {file_path}: {e}")
                    continue
        
        elapsed = time.time() - start_time
        print(f"Loaded {len(population_records)} samples in {elapsed:.2f}s")
        
        return {
            "success": True,
            "records": population_records
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/population/region/{region_str}")
async def get_population_region_data(region_str: str, folder_path: str):
    """Get population/cohort data for a specific region using parallel processing (legacy - loads all at once)"""
    try:
        folder = Path(folder_path)
        if not folder.exists():
            raise HTTPException(status_code=404, detail="Population folder not found")
        
        # Find all VCF files in the folder (.vcf.gz or .vcf)
        vcf_files = list(folder.glob("*.vcf.gz")) + list(folder.glob("*.vcf"))
        
        if not vcf_files:
            raise HTTPException(status_code=404, detail="No VCF files found in the folder")
        
        population_records = {}
        
        # Prepare arguments for parallel processing (3 args: file_path, region_str, cohort_mode)
        # Use None for cohort_mode to auto-detect format (individual mode has no explicit mode)
        file_args = [(str(vcf_file), region_str, None) for vcf_file in vcf_files]
        
        # Use ProcessPoolExecutor for true parallel processing (CPU-bound parsing)
        with ProcessPoolExecutor(max_workers=COHORT_WORKERS) as executor:
            # Submit all tasks
            future_to_file = {executor.submit(process_single_vcf_file, args): args[0] for args in file_args}
            
            # Collect results as they complete
            for future in as_completed(future_to_file):
                try:
                    result = future.result()
                    if result is not None:
                        sample_name, record = result
                        population_records[sample_name] = record
                except Exception as e:
                    file_path = future_to_file[future]
                    print(f"Error processing {file_path}: {e}")
                    continue
        
        # Return after processing all files
        return {
            "success": True,
            "records": population_records
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

# Pathogenic catalog cache
_pathogenic_catalog = None

def load_pathogenic_catalog():
    """Load pathogenic TR catalog from BED file - works with or without pandas"""
    global _pathogenic_catalog
    if _pathogenic_catalog is not None:
        return _pathogenic_catalog
    
    # try to find where the pathogenic catalog file is
    backend_dir = Path(__file__).parent
    catalog_paths = [
        backend_dir.parent.parent / "src" / "proletract" / "data" / "pathogenic_TRs.bed",
        backend_dir / "data" / "pathogenic_TRs.bed",
        Path("/confidential/home01/Calraei/tandemrepeats/ProleTRact/src/proletract/data/pathogenic_TRs.bed"),
    ]
    
    catalog_path = None
    for path in catalog_paths:
        if path.exists():
            catalog_path = path
            break
    
    if catalog_path is None:
        print("Warning: Pathogenic TR catalog not found. Tried paths:")
        for path in catalog_paths:
            print(f"  - {path} (exists: {path.exists()})")
        # try pandas if we have it, otherwise just use a list
        try:
            import pandas as pd
            _pathogenic_catalog = pd.DataFrame()
        except ImportError:
            _pathogenic_catalog = []
        return _pathogenic_catalog
    
    print(f"Loading pathogenic catalog from: {catalog_path}")
    
    # try pandas first since its faster
    try:
        import pandas as pd
        _pathogenic_catalog = pd.read_csv(catalog_path, sep="\t", header=None)
        _pathogenic_catalog.columns = ["chrom", "start", "end", "motif", "pathogenic_min", "inheritance", "disease", "gene"]
        _pathogenic_catalog["region"] = (
            _pathogenic_catalog["chrom"].astype(str) + ":" + 
            _pathogenic_catalog["start"].astype(str) + "-" + 
            _pathogenic_catalog["end"].astype(str)
        )
        print(f"Loaded pathogenic catalog with {len(_pathogenic_catalog)} regions (using pandas)")
        print(f"Sample regions: {_pathogenic_catalog[['chrom', 'start', 'end', 'gene']].head(3).to_string()}")
        return _pathogenic_catalog
    except ImportError:
        print("Pandas not available, loading catalog manually...")
    except Exception as e:
        print(f"Error loading with pandas: {e}, trying manual load...")
    
    # fallback: load it manually if pandas isnt available
    try:
        catalog_data = []
        with open(catalog_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= 8:
                    catalog_data.append({
                        'chrom': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'motif': parts[3],
                        'pathogenic_min': float(parts[4]) if parts[4] and parts[4] != '.' else None,
                        'inheritance': parts[5] if parts[5] else None,
                        'disease': parts[6] if parts[6] else None,
                        'gene': parts[7] if parts[7] else None,
                        'region': f"{parts[0]}:{parts[1]}-{parts[2]}"
                    })
        _pathogenic_catalog = catalog_data
        print(f"Loaded pathogenic catalog with {len(catalog_data)} regions (manual load)")
        if len(catalog_data) > 0:
            print(f"Sample region: {catalog_data[0]['region']} (gene: {catalog_data[0].get('gene', 'N/A')})")
        return _pathogenic_catalog
    except Exception as e:
        print(f"Error loading pathogenic catalog manually: {e}")
        _pathogenic_catalog = []
        return _pathogenic_catalog

@app.get("/api/pathogenic/check")
async def check_pathogenicity(chr: str, start: int, end: int):
    """Check if a region overlaps with pathogenic catalog"""
    try:
        print(f"Checking pathogenicity for {chr}:{start}-{end}")
        catalog = load_pathogenic_catalog()
        
        # handle both pandas dataframe and list of dicts
        if catalog is None:
            print(f"Catalog is None")
            return {
                "pathogenic": False
            }
        
        # check if its a list (means we loaded it manually without pandas)
        if isinstance(catalog, list):
            if len(catalog) == 0:
                print(f"Catalog is empty list")
                return {
                    "pathogenic": False
                }
            
            print(f"Catalog loaded: {len(catalog)} regions (manual load, no pandas)")
            
            # normalize chr name (add 'chr' prefix if missing)
            chr_normalized = chr if chr.startswith('chr') else f'chr{chr}'
            
            # filter by chromosome first
            chr_matches = [entry for entry in catalog 
                          if entry['chrom'] == chr_normalized or entry['chrom'] == chr]
            
            if len(chr_matches) == 0:
                print(f"No matches for chromosome {chr} (tried {chr_normalized} and {chr})")
                return {
                    "pathogenic": False
                }
            
            print(f"Found {len(chr_matches)} entries for chromosome {chr}")
            
            # find the best matching region (with some flexibility on coordinates)
            best_match = None
            best_distance = float('inf')
            
            for entry in chr_matches:
                start_diff = abs(entry['start'] - start)
                end_diff = abs(entry['end'] - end)
                total_distance = start_diff + end_diff
                
                # check if it matches (overlaps or is within 10bp)
                is_overlap = (entry['start'] <= end) and (entry['end'] >= start)
                is_close = (start_diff <= 10) and (end_diff <= 10)
                
                if (is_overlap or is_close) and total_distance < best_distance:
                    best_match = entry
                    best_distance = total_distance
            
            if best_match:
                pathogenic_min_val = best_match.get('pathogenic_min')
                pathogenic_threshold = int(float(pathogenic_min_val)) if pathogenic_min_val is not None else None
                
                print(f"Found match for {chr}:{start}-{end} with catalog entry {best_match['chrom']}:{best_match['start']}-{best_match['end']} (distance={best_distance}bp): threshold={pathogenic_threshold}, gene={best_match.get('gene', 'N/A')}")
                return {
                    "pathogenic": True,
                    "chr": str(best_match["chrom"]),
                    "start": int(best_match["start"]),
                    "end": int(best_match["end"]),
                    "gene": best_match.get("gene"),
                    "disease": best_match.get("disease"),
                    "inheritance": best_match.get("inheritance"),
                    "pathogenic_threshold": pathogenic_threshold,
                    "motif": best_match.get("motif")
                }
            
            print(f"No pathogenic match found for {chr}:{start}-{end} (checked {len(chr_matches)} entries)")
            return {
                "pathogenic": False
            }
        
        # Handle pandas DataFrame
        try:
            import pandas as pd
        except ImportError:
            print("Pandas not available and catalog is not a list")
            return {
                "pathogenic": False
            }
        
        if not isinstance(catalog, pd.DataFrame):
            print(f"Catalog is unexpected type: {type(catalog)}")
            return {
                "pathogenic": False
            }
        
        if catalog.empty:
            print(f"Catalog DataFrame is empty (shape: {catalog.shape})")
            return {
                "pathogenic": False
            }
        
        print(f"Catalog loaded: {len(catalog)} regions (pandas DataFrame)")
        
        # Normalize chromosome name (handle both 'chr13' and '13')
        chr_normalized = chr if chr.startswith('chr') else f'chr{chr}'
        
        # Check if region overlaps with any pathogenic entry
        # First try exact chromosome match
        chr_matches = catalog[catalog["chrom"] == chr_normalized]
        if len(chr_matches) == 0:
            # Try without 'chr' prefix
            chr_matches = catalog[catalog["chrom"] == chr]
        
        if len(chr_matches) == 0:
            print(f"No matches for chromosome {chr} (tried {chr_normalized} and {chr})")
            return {
                "pathogenic": False
            }
        
        print(f"Found {len(chr_matches)} entries for chromosome {chr}")
        
        # Check for overlap with tolerance for small coordinate differences
        # Make a copy to avoid SettingWithCopyWarning
        chr_matches = chr_matches.copy()
        chr_matches["start_diff"] = abs(chr_matches["start"] - start)
        chr_matches["end_diff"] = abs(chr_matches["end"] - end)
        chr_matches["total_distance"] = chr_matches["start_diff"] + chr_matches["end_diff"]
        
        # Find overlaps: standard overlap OR coordinates within 10bp
        overlaps = chr_matches[
            (
                # Standard overlap check
                (chr_matches["start"] <= end) & (chr_matches["end"] >= start)
            ) | (
                # Allow matches if start and end are both within 10bp
                (chr_matches["start_diff"] <= 10) & (chr_matches["end_diff"] <= 10)
            )
        ]
        
        if not overlaps.empty:
            # Find the best match (closest coordinates)
            best_match = overlaps.loc[overlaps["total_distance"].idxmin()]
            row = best_match
            
            # Handle pathogenic_min which might be float or None
            pathogenic_min_val = row["pathogenic_min"]
            if pd.notna(pathogenic_min_val):
                # Convert to int, handling float values
                pathogenic_threshold = int(float(pathogenic_min_val))
            else:
                pathogenic_threshold = None
            
            print(f"Found overlap for {chr}:{start}-{end} with catalog entry {row['chrom']}:{row['start']}-{row['end']} (distance={row['total_distance']}bp): threshold={pathogenic_threshold}, gene={row.get('gene', 'N/A')}")
            return {
                "pathogenic": True,
                "chr": str(row["chrom"]),
                "start": int(row["start"]),
                "end": int(row["end"]),
                "gene": str(row["gene"]) if pd.notna(row["gene"]) else None,
                "disease": str(row["disease"]) if pd.notna(row["disease"]) else None,
                "inheritance": str(row["inheritance"]) if pd.notna(row["inheritance"]) else None,
                "pathogenic_threshold": pathogenic_threshold,
                "motif": str(row["motif"]) if pd.notna(row["motif"]) else None
            }
        
        print(f"No pathogenic match found for {chr}:{start}-{end} (checked {len(chr_matches)} entries)")
        
        return {
            "pathogenic": False
        }
    except Exception as e:
        print(f"Error checking pathogenicity: {e}")
        return {
            "pathogenic": False
        }

@app.get("/api/pathogenic/search")
async def search_by_gene(gene: str):
    """Search for regions by gene name in pathogenic catalog"""
    try:
        print(f"Searching for gene: {gene}")
        catalog = load_pathogenic_catalog()
        
        if catalog is None:
            return {
                "success": False,
                "regions": [],
                "message": "Pathogenic catalog not available"
            }
        
        gene_lower = gene.lower().strip()
        matching_regions = []
        
        # Handle list format (manual load without pandas)
        if isinstance(catalog, list):
            if len(catalog) == 0:
                return {
                    "success": False,
                    "regions": [],
                    "message": "Pathogenic catalog is empty"
                }
            
            for entry in catalog:
                entry_gene = entry.get('gene', '')
                if entry_gene and gene_lower in entry_gene.lower():
                    matching_regions.append({
                        "region": entry.get('region', ''),
                        "chr": entry.get('chrom', ''),
                        "start": entry.get('start', 0),
                        "end": entry.get('end', 0),
                        "gene": entry_gene,
                        "disease": entry.get('disease'),
                        "inheritance": entry.get('inheritance'),
                        "motif": entry.get('motif'),
                        "pathogenic_threshold": entry.get('pathogenic_min')
                    })
        else:
            # Handle pandas DataFrame
            try:
                import pandas as pd
                if isinstance(catalog, pd.DataFrame) and len(catalog) > 0:
                    # Filter by gene (case-insensitive)
                    gene_matches = catalog[
                        catalog["gene"].astype(str).str.lower().str.contains(gene_lower, na=False)
                    ]
                    
                    for _, row in gene_matches.iterrows():
                        matching_regions.append({
                            "region": str(row.get("region", "")),
                            "chr": str(row.get("chrom", "")),
                            "start": int(row.get("start", 0)),
                            "end": int(row.get("end", 0)),
                            "gene": str(row.get("gene", "")) if pd.notna(row.get("gene")) else None,
                            "disease": str(row.get("disease", "")) if pd.notna(row.get("disease")) else None,
                            "inheritance": str(row.get("inheritance", "")) if pd.notna(row.get("inheritance")) else None,
                            "motif": str(row.get("motif", "")) if pd.notna(row.get("motif")) else None,
                            "pathogenic_threshold": float(row.get("pathogenic_min")) if pd.notna(row.get("pathogenic_min")) else None
                        })
            except ImportError:
                return {
                    "success": False,
                    "regions": [],
                    "message": "Pandas not available"
                }
        
        print(f"Found {len(matching_regions)} regions for gene '{gene}'")
        return {
            "success": True,
            "regions": matching_regions,
            "count": len(matching_regions)
        }
    except Exception as e:
        print(f"Error searching by gene: {e}")
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8502)
