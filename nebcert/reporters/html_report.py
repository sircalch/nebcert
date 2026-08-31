"""
Interactive HTML report dashboard generator for NEBCert.
"""

import os
import jinja2
from nebcert.core.scoring import ReactionPathwayReport

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NEBCert Reaction Kinetics & NEB Validation Report</title>
    <style>
        :root {
            --bg-color: #0f172a;
            --card-bg: #1e293b;
            --card-border: #334155;
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --pass-color: #10b981;
            --pass-bg: rgba(16, 185, 129, 0.15);
            --warn-color: #f59e0b;
            --warn-bg: rgba(245, 158, 11, 0.15);
            --fail-color: #ef4444;
            --fail-bg: rgba(239, 68, 68, 0.15);
            --accent-blue: #38bdf8;
            --accent-purple: #c084fc;
        }
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-primary);
            line-height: 1.6;
            padding: 2rem 1rem;
        }
        .container { max-width: 1200px; margin: 0 auto; }
        .header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            border-bottom: 1px solid var(--card-border);
            padding-bottom: 1.5rem;
            margin-bottom: 2rem;
            flex-wrap: wrap;
            gap: 1rem;
        }
        .title-group h1 { font-size: 2rem; font-weight: 700; }
        .title-group p { color: var(--text-secondary); font-size: 0.95rem; }
        .status-badge {
            display: inline-flex;
            align-items: center;
            padding: 0.5rem 1.25rem;
            border-radius: 9999px;
            font-weight: 700;
            font-size: 1.1rem;
            letter-spacing: 0.05em;
            text-transform: uppercase;
        }
        .badge-pass { background-color: var(--pass-bg); color: var(--pass-color); border: 1px solid var(--pass-color); }
        .badge-warning { background-color: var(--warn-bg); color: var(--warn-color); border: 1px solid var(--warn-color); }
        .badge-fail { background-color: var(--fail-bg); color: var(--fail-color); border: 1px solid var(--fail-color); }

        .grid-cards {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 1.25rem;
            margin-bottom: 2rem;
        }
        .card {
            background-color: var(--card-bg);
            border: 1px solid var(--card-border);
            border-radius: 0.75rem;
            padding: 1.25rem;
            box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.2);
        }
        .card-label { font-size: 0.8rem; text-transform: uppercase; color: var(--text-secondary); margin-bottom: 0.4rem; }
        .card-value { font-size: 1.4rem; font-weight: 700; color: var(--text-primary); }
        .card-subtext { font-size: 0.85rem; color: var(--text-secondary); margin-top: 0.25rem; }

        .section-title { font-size: 1.3rem; font-weight: 600; margin-bottom: 1rem; color: var(--accent-blue); }
        table {
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 2rem;
            font-size: 0.95rem;
            background-color: var(--card-bg);
            border-radius: 0.75rem;
            overflow: hidden;
            border: 1px solid var(--card-border);
        }
        th, td { padding: 0.85rem 1rem; text-align: left; border-bottom: 1px solid var(--card-border); }
        th { background-color: rgba(255, 255, 255, 0.03); color: var(--text-secondary); font-weight: 600; text-transform: uppercase; font-size: 0.85rem; }
        tr:hover td { background-color: rgba(255, 255, 255, 0.02); }

        .tag { display: inline-block; padding: 0.2rem 0.6rem; border-radius: 0.375rem; font-size: 0.75rem; font-weight: 700; }
        .tag-pass { background-color: var(--pass-bg); color: var(--pass-color); }
        .tag-warning { background-color: var(--warn-bg); color: var(--warn-color); }
        .tag-fail { background-color: var(--fail-bg); color: var(--fail-color); }

        .box { background-color: var(--card-bg); border: 1px solid var(--card-border); border-radius: 0.75rem; padding: 1.25rem; margin-bottom: 2rem; }
        pre { background-color: rgba(0, 0, 0, 0.4); padding: 1rem; border-radius: 0.5rem; color: #38bdf8; font-family: monospace; font-size: 0.85rem; white-space: pre-wrap; }
        .btn-copy { background-color: #2563eb; color: white; border: none; padding: 0.4rem 0.8rem; border-radius: 0.375rem; cursor: pointer; font-size: 0.8rem; margin-top: 0.5rem; }
        .btn-copy:hover { background-color: #1d4ed8; }

        footer { text-align: center; font-size: 0.85rem; color: var(--text-secondary); margin-top: 3rem; padding-top: 1.5rem; border-top: 1px solid var(--card-border); }
    </style>
</head>
<body>
    <div class="container">
        <header class="header">
            <div class="title-group">
                <h1>NEBCert Reaction Kinetics & NEB Report</h1>
                <p>{{ report.metadata.reaction }} &bull; {{ report.metadata.functional }} ({{ report.metadata.software }})</p>
            </div>
            <div>
                <span class="status-badge badge-{{ report.overall_status.lower() }}">
                    {{ report.overall_status }}
                </span>
            </div>
        </header>

        <div class="grid-cards">
            <div class="card">
                <div class="card-label">Activation Barrier (E_a^fwd)</div>
                <div class="card-value">
                    {% if report.neb_profile %}
                    {{ "%.2f"|format(report.neb_profile.e_forward_barrier_kcal_mol) }} kcal/mol
                    {% else %}
                    N/A
                    {% endif %}
                </div>
                <div class="card-subtext">
                    {% if report.neb_profile %}
                    {{ "%.3f"|format(report.neb_profile.e_forward_barrier_ev) }} eV (&Delta;E_rxn = {{ "%.2f"|format(report.neb_profile.delta_e_reaction_kcal_mol) }} kcal/mol)
                    {% else %}
                    No NEB profile provided
                    {% endif %}
                </div>
            </div>
            <div class="card">
                <div class="card-label">TS Imaginary Frequency</div>
                <div class="card-value">
                    {% if report.ts_frequency and report.ts_frequency.imaginary_frequency_cm1 %}
                    {{ "%.1f"|format(report.ts_frequency.imaginary_frequency_cm1) }} cm&sup1;
                    {% else %}
                    N/A
                    {% endif %}
                </div>
                <div class="card-subtext">
                    {% if report.ts_frequency %}
                    {{ report.ts_frequency.n_imaginary_frequencies }} imaginary mode(s)
                    {% else %}
                    No frequency data
                    {% endif %}
                </div>
            </div>
            <div class="card">
                <div class="card-label">TST Rate Constant (298 K)</div>
                <div class="card-value">
                    {% if report.tst_kinetics %}
                    {{ "%.2e"|format(report.tst_kinetics.k_298_s_minus_1) }} s&sup1;
                    {% else %}
                    N/A
                    {% endif %}
                </div>
                <div class="card-subtext">
                    {% if report.tst_kinetics %}
                    Half-life: {{ "%.2e"|format(report.tst_kinetics.half_life_298_s) }} s
                    {% else %}
                    No kinetics computed
                    {% endif %}
                </div>
            </div>
            <div class="card">
                <div class="card-label">Quantum Tunneling (&kappa;_298)</div>
                <div class="card-value">
                    {% if report.tunneling %}
                    {{ "%.2f"|format(report.tunneling.kappa_eckart_298) }}
                    {% else %}
                    N/A
                    {% endif %}
                </div>
                <div class="card-subtext">
                    {% if report.tunneling %}
                    Wigner factor: {{ "%.2f"|format(report.tunneling.kappa_wigner_298) }}
                    {% else %}
                    No tunneling computed
                    {% endif %}
                </div>
            </div>
        </div>

        <h2 class="section-title">Reaction Pathway & Transition State Verification Criteria</h2>
        <table>
            <thead>
                <tr>
                    <th>Reaction Property / Criterion</th>
                    <th>Measured Value</th>
                    <th>Standard Recommended Threshold</th>
                    <th>Status</th>
                    <th>Diagnostic Interpretation</th>
                </tr>
            </thead>
            <tbody>
                {% if report.neb_profile %}
                <tr>
                    <td><strong>Minimum Energy Path (MEP)</strong></td>
                    <td>E_a^fwd = {{ "%.2f"|format(report.neb_profile.e_forward_barrier_kcal_mol) }} kcal/mol ({{ report.neb_profile.n_images }} images)</td>
                    <td>Continuous single-step MEP</td>
                    <td><span class="tag tag-{{ report.neb_profile.status.lower() }}">{{ report.neb_profile.status }}</span></td>
                    <td style="color: var(--text-secondary)">{{ report.neb_profile.diagnostic_message }}</td>
                </tr>
                {% endif %}
                {% if report.ts_frequency %}
                <tr>
                    <td><strong>Transition State 1st-Order Saddle Point</strong></td>
                    <td>{{ report.ts_frequency.n_imaginary_frequencies }} imaginary mode ({{ "%.1f"|format(report.ts_frequency.imaginary_frequency_cm1 if report.ts_frequency.imaginary_frequency_cm1 else 0.0) }} cm&sup1;)</td>
                    <td>Strictly 1 imaginary frequency (|&nu;| &ge; 50 cm&sup1;)</td>
                    <td><span class="tag tag-{{ report.ts_frequency.status.lower() }}">{{ report.ts_frequency.status }}</span></td>
                    <td style="color: var(--text-secondary)">{{ report.ts_frequency.diagnostic_message }}</td>
                </tr>
                {% endif %}
                {% if report.tst_kinetics %}
                <tr>
                    <td><strong>Eyring TST Rate Constants</strong></td>
                    <td>k(298 K) = {{ "%.2e"|format(report.tst_kinetics.k_298_s_minus_1) }} s&sup1; (E_a = {{ "%.2f"|format(report.tst_kinetics.arrhenius_e_activation_kcal_mol) }} kcal/mol)</td>
                    <td>Eyring-Polanyi TST</td>
                    <td><span class="tag tag-{{ report.tst_kinetics.status.lower() }}">{{ report.tst_kinetics.status }}</span></td>
                    <td style="color: var(--text-secondary)">{{ report.tst_kinetics.diagnostic_message }}</td>
                </tr>
                {% endif %}
                {% if report.tunneling %}
                <tr>
                    <td><strong>Eckart Quantum Tunneling</strong></td>
                    <td>&kappa;_Eckart = {{ "%.2f"|format(report.tunneling.kappa_eckart_298) }} (298 K)</td>
                    <td>1D Asymmetric Eckart Transmission</td>
                    <td><span class="tag tag-{{ report.tunneling.status.lower() }}">{{ report.tunneling.status }}</span></td>
                    <td style="color: var(--text-secondary)">{{ report.tunneling.diagnostic_message }}</td>
                </tr>
                {% endif %}
            </tbody>
        </table>

        {% if report.recommendations %}
        <div class="box" style="border-left: 4px solid var(--warn-color);">
            <h3 style="color: var(--warn-color); margin-bottom: 0.5rem;">Diagnostic Recommendations & Methodological Alerts</h3>
            <ul style="padding-left: 1.25rem;">
                {% for rec in report.recommendations %}
                <li style="margin-bottom: 0.25rem; color: var(--text-secondary);">{{ rec }}</li>
                {% endfor %}
            </ul>
        </div>
        {% endif %}

        <h2 class="section-title">Computational Methods (Publication Ready)</h2>
        <div class="box">
            <pre id="methodsSnippet">{{ methods_text }}</pre>
            <button class="btn-copy" onclick="copyToClipboard('methodsSnippet')">Copy Methods Snippet</button>
        </div>

        <h2 class="section-title">BibTeX Citation</h2>
        <div class="box">
            <pre id="bibSnippet">{{ citation_bib }}</pre>
            <button class="btn-copy" onclick="copyToClipboard('bibSnippet')">Copy BibTeX</button>
        </div>

        <footer>
            Generated automatically by <strong>NEBCert v1.0.0</strong> &bull; Reaction Kinetics & Transition State Certification &bull; Monreal-Hernández, 2026.
        </footer>
    </div>

    <script>
        function copyToClipboard(elementId) {
            const text = document.getElementById(elementId).innerText;
            navigator.clipboard.writeText(text).then(() => {
                alert('Copied to clipboard!');
            }).catch(err => {
                console.error('Error copying: ', err);
            });
        }
    </script>
</body>
</html>
"""


def generate_nebcert_html_report(
    report: ReactionPathwayReport,
    output_path: str,
    methods_text: str = "",
    citation_bib: str = ""
) -> str:
    """
    Renders HTML report template and writes it to disk.
    """
    template = jinja2.Template(HTML_TEMPLATE)
    rendered = template.render(
        report=report,
        methods_text=methods_text,
        citation_bib=citation_bib
    )
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(rendered)
    return output_path
