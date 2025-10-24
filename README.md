<h1>🧬 ProleTRact</h1>
<p>This repository contains a <strong>Tandem Repeat Visualization Tool</strong> designed using Streamlit. The tool processes Variant Call Format (VCF) files, visualizing tandem repeats in an intuitive, interactive format. Users can explore motifs, compare alleles, and gain insights into the structure of tandem repeats, enhancing their ability to interpret genomic variation.</p>

<h2>Features</h2>
<ul>
  <li><strong>Chohort and indivisual mode</li>
  <li><strong>Dynamic Sequence Visualization:</strong> Displays the sequence with color-coded motifs, highlighting interruptions in red and making it easy to identify motif patterns and interruptions.</li>
  <li><strong>Motif Comparison Across Alleles:</strong> Visualizes motifs for each allele separately, providing a side-by-side comparison of motif structures and copy numbers.</li>
  <li><strong>Population comparsion</li>

</ul>

<h2>Installation</h2>
<ol>
  <li>Clone this repository:
    <pre><code>git clone git@github.com:Lionward/tandemtwister-vis.git</code></pre>
  </li>
  <li>Navigate to the project directory:
    <pre><code>cd tandemtwister-vis</code></pre>
  </li>
  <li>Install the required dependencies:
  </li>
</ol>

<h2>Usage</h2>
<ol>
  <li>Run the Streamlit app:</li>
    <pre><code>streamlit run app.py</code></pre>
    <pre><code> ssh -L 8501:localhost:8501 phaedra</code></pre>
  </li>
  <li>Upload your VCF file using the file uploader in the sidebar.</li>
  <li>Navigate through different tandem repeat regions using the "Next" and "Previous" buttons or input the region manually (e.g., <code>chr1:1000-2000</code>).</li>
</ol>

<h2>Dependencies</h2>
<ul>
  <li>Streamlit</li>
  <li>pysam</li>
  <li>matplotlib</li>
  <li>pandas</li>
  <li>altair</li>
  <li>plotly</li>
</ul>

<h2></h2>
<p>Below is an example screenshot of the Tandem Repeat Visualization Tool in action:</p>
<img src="EXAMPLE.png" alt="Tandem Repeat Visualization Example" style="max-width: 100%; height: auto; border: 1px solid #ccc; padding: 10px;">

<h2>Contributing</h2>
<p>If you'd like to contribute to this project, please fork the repository, make your changes, and submit a pull request. Any contributions are welcome!</p>

<h2>License</h2>
<p>This project is licensed under the MIT License - see the <code>LICENSE</code> file for details.</p>
