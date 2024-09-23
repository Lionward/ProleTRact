<h1>Tandem Repeat Visualization Tool</h1>

<p>This repository provides a <strong>Tandem Repeat Visualization Tool</strong> developed using Streamlit. The tool processes variant call files (VCF) that are produced by tandemtiwster to visualize tandem repeats in a user-friendly and interactive way. It allows users to explore motifs in a specific genomic region, offering zoom-in functionality, motif highlighting, and summary statistics.</p>

<h2>Features</h2>
<ul>
  <li><strong>Motif Visualization:</strong> Visualizes motif sequences across alleles with distinct colors and labeled interruptions.</li>
  <li><strong>Allele Comparison:</strong> Displays motif occurrences for each allele separately and compares them side by side.</li>
  <li><strong>Interactive Zooming:</strong> Allows users to zoom in on specific regions of the sequence.</li>
  <li><strong>Hover Functionality:</strong> Highlights motifs when hovering over the legend, linking motifs in the sequence.</li>
  <li><strong>Summary Statistics:</strong> Displays total motif counts and records in the sidebar.</li>
</ul>

<h2>Installation</h2>
<ol>
  <li>Clone this repository:</li>
  <pre><code>git clone git@github.com:Lionward/tandemtwister-vis.git</code></pre>
  <li>Navigate to the project directory:</li>
  <pre><code>cd tandemtwister-vis</code></pre>
</ol>

<h2>Usage</h2>
<ol>
  <li>Run the Streamlit app:</li>
  <pre><code>streamlit run visualize_TRs.py</code></pre>
  <pre><code> ssh -L 8501:localhost:8501 phaedra</code></pre>
  <li>Upload your VCF file using the file uploader in the sidebar.</li>
  <li>Use the zoom slider to adjust the visualized region of the sequence.</li>
  <li>Navigate through different tandem repeat regions using the "Next" and "Previous" buttons or you can manually input the region e.g. chr1:1000-2000.</li>
</ol>


<h2>Dependencies</h2>
<ul>
  <li>Streamlit</li>
  <li>pysam</li>
  <li>plotly</li>
  <li>matplotlib</li>
  <li>pandas</li>
</ul>

<h2>Contributing</h2>
<p>If you'd like to contribute to this project, please fork the repository, make your changes, and submit a pull request. Any contributions are welcome!</p>

<h2>License</h2>
<p>This project is licensed under the MIT License - see the <code>LICENSE</code> file for details.</p>
