export interface VCFData {
  path: string;
  totalRegions: number | null;
  availableGenotypes: string[];
  availableChromosomes?: string[];
}

export interface RegionInfo {
  id: string;
  region: string;
  genotype: string;
}

export interface FilterResponse {
  records: RegionInfo[];
  total_matching: number;
  total_regions: number;
  current_page: number;
  total_pages: number;
  available_genotypes?: string[];
}

