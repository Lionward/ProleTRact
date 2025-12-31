import React, { useMemo } from 'react';
import './ExternalLinksPanel.css';
import externalLinksData from '../data/externalLinks.json';

interface ExternalLinksPanelProps {
  chr?: string;
  pos?: number;
  stop?: number;
  gene?: string;
  region?: string;
}

const ExternalLinksPanel: React.FC<ExternalLinksPanelProps> = ({
  chr,
  pos,
  stop,
  gene,
  region
}) => {
  const formatUrl = (urlTemplate: string): string => {
    if (!chr || !pos || !stop) return urlTemplate;
    
    // Normalize chromosome name (remove 'chr' prefix if present, then add it back for URLs)
    const chrNormalized = chr.startsWith('chr') ? chr : `chr${chr}`;
    const chrForUrl = chrNormalized.replace('chr', '');
    
    return urlTemplate
      .replace(/{chr}/g, chrForUrl)
      .replace(/{pos}/g, pos.toString())
      .replace(/{stop}/g, stop.toString())
      .replace(/{gene}/g, gene || '');
  };

  const links = useMemo(() => {
    const allLinks: Array<{
      category: string;
      name: string;
      url: string;
      icon: string;
      description: string;
    }> = [];

    // Genome Browsers
    externalLinksData.genomeBrowsers.forEach(link => {
      if (chr && pos && stop) {
        allLinks.push({
          category: 'Genome Browsers',
          ...link,
          url: formatUrl(link.url)
        });
      }
    });

    // Variant Databases
    externalLinksData.variantDatabases.forEach(link => {
      if (chr && pos) {
        allLinks.push({
          category: 'Variant Databases',
          ...link,
          url: formatUrl(link.url)
        });
      }
    });

    // Tandem Repeat Databases (no region-specific URL needed)
    externalLinksData.tandemRepeatDatabases.forEach(link => {
      allLinks.push({
        category: 'Tandem Repeat Databases',
        ...link,
        url: link.url
      });
    });

    // Population Databases (no region-specific URL needed)
    externalLinksData.populationDatabases.forEach(link => {
      allLinks.push({
        category: 'Population Databases',
        ...link,
        url: link.url
      });
    });

    // Gene Annotation
    if (gene) {
      externalLinksData.geneAnnotation.forEach(link => {
        allLinks.push({
          category: 'Gene Annotation',
          ...link,
          url: formatUrl(link.url)
        });
      });
    }

    return allLinks;
  }, [chr, pos, stop, gene]);

  const groupedLinks = useMemo(() => {
    const groups: { [key: string]: typeof links } = {};
    links.forEach(link => {
      if (!groups[link.category]) {
        groups[link.category] = [];
      }
      groups[link.category].push(link);
    });
    return groups;
  }, [links]);

  if (links.length === 0) {
    return null;
  }

  return (
    <div className="external-links-panel">
      <div className="external-links-header">
        <h3>🔗 External Resources</h3>
        <p className="external-links-subtitle">Quick links to related databases and tools</p>
      </div>
      
      <div className="external-links-content">
        {Object.entries(groupedLinks).map(([category, categoryLinks]) => (
          <div key={category} className="external-links-category">
            <h4 className="external-links-category-title">{category}</h4>
            <div className="external-links-grid">
              {categoryLinks.map((link, index) => (
                <a
                  key={`${category}-${index}`}
                  href={link.url}
                  target="_blank"
                  rel="noopener noreferrer"
                  className="external-link-item"
                  title={link.description}
                >
                  <span className="external-link-icon">{link.icon}</span>
                  <span className="external-link-name">{link.name}</span>
                  <span className="external-link-arrow">↗</span>
                </a>
              ))}
            </div>
          </div>
        ))}
      </div>
      
      {region && (
        <div className="external-links-region-info">
          <small>Region: {region}</small>
        </div>
      )}
    </div>
  );
};

export default ExternalLinksPanel;




