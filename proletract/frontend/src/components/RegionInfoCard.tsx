import React, { useMemo } from 'react';
import './PopulationComparison.css'; // Reuse styles

interface RegionInfoCardProps {
  chr: string;
  pos: number;
  stop: number;
  region?: string;
  currentAlleles?: string[];
}

const RegionInfoCard: React.FC<RegionInfoCardProps> = ({
  chr,
  pos,
  stop,
  region,
  currentAlleles = []
}) => {
  // Generate database URLs for the current region
  const databaseUrls = useMemo(() => {
    const urls: { [key: string]: { url: string; icon: string; color: string } } = {};
    
    const start = pos;
    const end = stop;
    const pos_range = `${start}-${end}`;
    
    // Remove 'chr' prefix if present for some databases
    const chrom_no_chr = chr.replace(/^chr/i, '');
    
    if (chr && pos_range) {
      // UCSC
      urls['UCSC'] = {
        url: `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${chr}:${pos_range}`,
        icon: '🌐',
        color: 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)'
      };
      
      if (chrom_no_chr) {
        // Ensembl
        urls['Ensembl'] = {
          url: `https://www.ensembl.org/Homo_sapiens/Location/View?r=${chrom_no_chr}:${pos_range}`,
          icon: '🧬',
          color: 'linear-gradient(135deg, #FF6B6B 0%, #EE5A6F 100%)'
        };
        
        // NCBI
        urls['NCBI'] = {
          url: start && end 
            ? `https://www.ncbi.nlm.nih.gov/genome/gdv/browser/?context=genome&acc=GCF_000001405.40&chr=${chrom_no_chr}&from=${start}&to=${end}`
            : `https://www.ncbi.nlm.nih.gov/genome/gdv/browser/?context=genome&acc=GCF_000001405.40&chr=${chrom_no_chr}`,
          icon: '📊',
          color: 'linear-gradient(135deg, #96CEB4 0%, #6C9A8B 100%)'
        };
        
        // DECIPHER
        urls['DECIPHER'] = {
          url: `https://www.deciphergenomics.org/browser#q/grch37:${chrom_no_chr}:${pos_range}`,
          icon: '🔍',
          color: 'linear-gradient(135deg, #A8E6CF 0%, #7FCDCD 100%)'
        };
      }
      
      if (chrom_no_chr && start && end) {
        // gnomAD
        urls['gnomAD'] = {
          url: `https://gnomad.broadinstitute.org/region/${chrom_no_chr}-${start}-${end}`,
          icon: '📈',
          color: 'linear-gradient(135deg, #FFE66D 0%, #FF6B6B 100%)'
        };
      }
      
      if (chr && pos_range) {
        // TRExplorer
        urls['TRExplorer'] = {
          url: `https://trexplorer.broadinstitute.org/#sc=isPathogenic&sd=DESC&showRs=1&searchQuery=${chr}:${pos_range}&showColumns=0i1i2i3i4i7i21i17`,
          icon: '🔬',
          color: 'linear-gradient(135deg, #10B981 0%, #059669 100%)'
        };
      }
    }
    
    return urls;
  }, [chr, pos, stop]);

  // Format region as chr:start-end (without commas)
  const regionString = `${chr}:${pos}-${stop}`;
  const length = stop - pos + 1;

  return (
    <div className="region-info-card-centered">
      <div className="region-display-center">
        <div className="region-main-display">
          <div className="region-string">{regionString}</div>
          <div className="region-length">{length.toLocaleString()} bp</div>
        </div>
        {currentAlleles.length > 0 && (
          <div className="region-alleles-center">
            {currentAlleles.map((allele, idx) => (
              <span key={idx} className="allele-badge-center">{allele}</span>
            ))}
          </div>
        )}
      </div>
      
      {/* Database Links Tags */}
      {Object.keys(databaseUrls).length > 0 && (
        <div className="database-links-container">
          <div className="database-links-label">View in External Databases:</div>
          <div className="database-links-tags">
            {Object.entries(databaseUrls).map(([name, { url, icon, color }]) => (
              <a
                key={name}
                href={url}
                target="_blank"
                rel="noopener noreferrer"
                className="database-link-tag"
                style={{ background: color }}
              >
                <span className="database-link-icon">{icon}</span>
                <span className="database-link-name">{name}</span>
              </a>
            ))}
          </div>
        </div>
      )}
    </div>
  );
};

export default RegionInfoCard;

