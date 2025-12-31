import React from 'react';
import './FeaturesPage.css';

const FeaturesPage: React.FC = () => {
  return (
    <div className="features-page">
      <div className="features-header">
        <h1>🚀 Upcoming Features & Roadmap</h1>
        <p className="features-subtitle">
          Discover what's coming next to ProleTRact
        </p>
      </div>

      <div className="features-content">
        {/* TandemTwister Link Section */}
        <section className="tandtwister-section">
          <div className="tandtwister-card">
            <div className="logos-container">
              <div className="logo-item">
                <img 
                  src="/ProleTRact_logo.svg" 
                  alt="ProleTRact Logo" 
                  className="companion-logo"
                />
              </div>
              <div className="logo-divider">+</div>
              <div className="logo-item">
                <img 
                  src="https://raw.githubusercontent.com/Lionward/tandemtwister/main/logo_tandemtwister.svg" 
                  alt="TandemTwister Logo" 
                  className="companion-logo"
                />
              </div>
            </div>
            <p className="companion-description">
              ProleTRact is the companion visualization tool for TandemTwister.
            </p>
            <a 
              href="https://github.com/Lionward/tandemtwister" 
              target="_blank" 
              rel="noopener noreferrer"
              className="tandtwister-link"
            >
              <span>🔗</span>
              <span>Visit TandemTwister on GitHub</span>
              <span>→</span>
            </a>
          </div>
        </section>

        {/* Upcoming Features Section */}
        <section className="features-section">
          <h2>📋 Upcoming Features</h2>
          
          <div className="feature-category">
            <h3>🔬 Enhanced Visualizations</h3>
            <div className="feature-list">
              <div className="feature-item">
                <div className="feature-status planned">Planned</div>
                <div className="feature-content">
                  <h4>IGV-like Read Browser</h4>
                  <p>Interactive read alignment viewer with individual read tracks, depth visualization, 
                  and zoom/pan functionality for detailed read-level analysis.</p>
                </div>
              </div>
            </div>
          </div>

          <div className="feature-category">
            <h3>📊 Advanced Analysis</h3>
            <div className="feature-list">
              <div className="feature-item">
                <div className="feature-status planned">Planned</div>
                <div className="feature-content">
                  <h4>HGSVC/Population Comparison</h4>
                  <p>Compare sample alleles against HGSVC population data to identify rare or novel 
                  expansions with allele frequency information.</p>
                </div>
              </div>
              <div className="feature-item">
                <div className="feature-status planned">Planned</div>
                <div className="feature-content">
                  <h4>Advanced Clustering</h4>
                  <p>Enhanced clustering algorithms with multiple distance metrics, silhouette score 
                  analysis, and interactive cluster exploration.</p>
                </div>
              </div>
            </div>
          </div>

          <div className="feature-category">
            <h3>⚡ Performance & Usability</h3>
            <div className="feature-list">
              <div className="feature-item">
                <div className="feature-status planned">Planned</div>
                <div className="feature-content">
                  <h4>Improved Large Cohort Handling</h4>
                  <p>Optimized data loading and rendering for large cohort datasets with virtual scrolling 
                  and progressive data loading.</p>
                </div>
              </div>
              <div className="feature-item">
                <div className="feature-status planned">Planned</div>
                <div className="feature-content">
                  <h4>Export Enhancements</h4>
                  <p>Additional export formats including PDF reports, batch export capabilities, and 
                  customizable export templates.</p>
                </div>
              </div>
              <div className="feature-item">
                <div className="feature-status planned">Planned</div>
                <div className="feature-content">
                  <h4>Session Management</h4>
                  <p>Save and restore analysis sessions with bookmarks, filters, and visualization states.</p>
                </div>
              </div>
            </div>
          </div>

          <div className="feature-category">
            <h3>🔗 Integration & API</h3>
            <div className="feature-list">
              <div className="feature-item">
                <div className="feature-status planned">Planned</div>
                <div className="feature-content">
                  <h4>REST API</h4>
                  <p>Programmatic access to ProleTRact functionality via REST API for integration with 
                  other bioinformatics tools and workflows.</p>
                </div>
              </div>
              <div className="feature-item">
                <div className="feature-status planned">Planned</div>
                <div className="feature-content">
                  <h4>Database Integration</h4>
                  <p>Direct integration with genomic databases for automatic annotation and metadata retrieval.</p>
                </div>
              </div>
            </div>
          </div>
        </section>

        {/* Current Features Section */}
        <section className="features-section">
          <h2>✅ Currently Available Features</h2>
          <div className="current-features-grid">
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Statistics Dashboard</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Region Visualization</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Population Comparison</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Cohort Analysis</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Genotype Filtering</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Stack Plot & Heatmap</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Bar Plot for Motif Count</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Cluster Plot</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Data Export (CSV, JSON, FASTA)</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Bookmarks</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Keyboard Shortcuts</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Theme Toggle (Dark/Light)</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Pathogenic TR Catalog Integration</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Pathogenic Threshold Visualization</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Pathogenic Region Filter</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Gene Search</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Pathogenicity Panel</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Allele Color Coding (Pathogenic/Non-pathogenic)</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Sequence Display Modes (Sequence & Bar View)</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>Scroll Position Preservation</span>
            </div>
            <div className="current-feature-item">
              <span className="feature-icon">✓</span>
              <span>View Mode Persistence</span>
            </div>
          </div>
        </section>

        {/* Contribution Section */}
        <section className="features-section">
          <h2>🤝 Contributing</h2>
          <p className="contribution-text">
            ProleTRact is an open-source project. Contributions, suggestions, and feedback are welcome! 
            Visit our GitHub repository to report issues, suggest features, or contribute code.
          </p>
        </section>
      </div>
    </div>
  );
};

export default FeaturesPage;

