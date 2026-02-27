import streamlit as st
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio import SeqIO, Entrez
from Bio.SeqUtils import MeltingTemp as mt
import matplotlib.pyplot as plt
import io, ssl, os

# 1. 보안 및 환경 설정
ssl._create_default_https_context = ssl._create_unverified_context
os.environ['CURL_CA_BUNDLE'] = ''
Entrez.email = "your_email@example.com" # NCBI 에티켓용

# 페이지 설정
st.set_page_config(page_title="Bio-Research Station", layout="wide")

# 커스텀 CSS (UI 디자인)
st.markdown("""
    <style>
    .stApp { background-color: #ffffff; }
    .main-title { font-size: 30px; font-weight: 800; color: #2D3436; margin-bottom: 25px; border-bottom: 2px solid #F1F2F6; padding-bottom: 10px; }
    .result-card { background-color: #FDFDFD; border: 1px solid #EAEAEA; padding: 20px; border-radius: 8px; margin-bottom: 15px; box-shadow: 0 2px 4px rgba(0,0,0,0.02); }
    .target-label { background-color: #2D3436; color: white; padding: 4px 12px; border-radius: 4px; font-weight: 600; font-size: 12px; display: inline-block; margin-bottom: 10px; }
    .guide-box { background-color: #F8F9FA; border-left: 4px solid #D1D4D9; padding: 15px; color: #636E72; font-size: 14px; }
    </style>
    """, unsafe_allow_html=True)

# 2. 로그인 로직
if 'auth' not in st.session_state:
    st.session_state['auth'] = False

def check_password():
    if st.session_state["pw_input"] == "knu2026":
        st.session_state['auth'] = True
        del st.session_state["pw_input"]
    else:
        st.error("비밀번호가 틀렸습니다.")

if not st.session_state['auth']:
    st.markdown('<div style="text-align:center; margin-top:150px;"><h1> Restricted Access</h1><p>허가된 사용자 전용 시스템입니다.</p></div>', unsafe_allow_html=True)
    st.text_input("Access Password", type="password", on_change=check_password, key="pw_input")
    st.stop()

# 3. 유전자 기능 분석 사전 정의
def get_pwn_analysis(title):
    t = title.lower()
    if any(w in t for w in ["cellulase", "pectate", "eng-1", "pel1", "expansin"]): return "식물 세포벽 분해 효소: 조직 침입 관여"
    if any(w in t for w in ["v-atpase", "hsp90", "ribosomal", "atp synthase"]): return "핵심 생존 유전자: 생명 유지 (RNAi 타겟)"
    if any(w in t for w in ["cytochrome", "cyp-", "gst-"]): return "약제 저항성: 살선충제 방어 기제"
    return "기능 분석 필요: 상세 정보 및 논문을 확인하세요."

# 4. 메인 레이아웃
st.markdown('<p class="main-title">🌲 B. xylophilus Integrated Research Station</p>', unsafe_allow_html=True)

tab1, tab2, tab3, tab4, tab5 = st.tabs(["🌲 Target Analysis", "🧬 Central Dogma", "📂 Converter", "🔍 Primer Design", "📊 Visualization"])

# --- 탭 1: 소나무재선충 분석 및 논문 검색 ---
with tab1:
    st.subheader("Target Search & Literature Analysis")
    col_in, col_guide = st.columns([3, 1])
    with col_in:
        seq_in = st.text_area("DNA Sequence (Target Analysis)", height=150, placeholder="분석할 서열을 입력하세요.")
        if st.button("RUN FULL ANALYSIS", use_container_width=True):
            if len(seq_in) < 15: st.warning("서열이 너무 짧습니다.")
            else:
                with st.spinner("NCBI BLAST & PubMed 검색 중..."):
                    try:
                        res = NCBIWWW.qblast("blastn", "nt", seq_in, entrez_query="Bursaphelenchus xylophilus [ORGN]")
                        rec = NCBIXML.read(res)
                        if rec.alignments:
                            for i, aln in enumerate(rec.alignments[:3]):
                                st.markdown(f"""
                                <div class="result-card">
                                    <div class="target-label">Candidate {i+1}</div>
                                    <div style="font-size:18px; font-weight:700;">
                                        <a href="https://www.ncbi.nlm.nih.gov/nuccore/{aln.accession}" target="_blank" style="text-decoration:none; color:#0984e3;">
                                            🔍 {aln.accession} (NCBI Link)
                                        </a>
                                    </div>
                                    <p style="font-size:13px; color:#636E72;">{aln.title[:150]}...</p>
                                    <p><b>기능 분석 결과:</b> {get_pwn_analysis(aln.title)}</p>
                                </div>""", unsafe_allow_html=True)
                                # 논문 검색
                                try:
                                    h = Entrez.esearch(db="pubmed", term=f"{aln.accession} Bursaphelenchus", retmax=2)
                                    pids = Entrez.read(h)["IdList"]
                                    if pids:
                                        st.write("📚 **Related Publications:**")
                                        for pid in pids:
                                            smh = Entrez.esummary(db="pubmed", id=pid); sm = Entrez.read(smh)
                                            st.markdown(f"- [{sm[0]['Title']}](https://pubmed.ncbi.nlm.nih.gov/{pid}/)")
                                except: pass
                                st.markdown("---")
                        else: st.info("결과가 없습니다.")
                    except Exception as e: st.error(f"오류 발생: {e}")
    with col_guide:
        st.markdown('<div class="guide-box"><b>BLAST Guide</b><br>재선충 유전자를 분석하고 관련 논문과 NCBI 원문 링크를 자동으로 매칭합니다.</div>', unsafe_allow_html=True)

# --- 탭 2: 전사/번역 ---
with tab2:
    st.subheader("Transcription & Translation")
    dna_val = st.text_area("DNA Sequence", key="dogma_in").upper().strip()
    if dna_val:
        s = Seq(dna_val)
        c1, c2 = st.columns(2)
        with c1: st.write("**RNA Sequence:**"); st.code(s.transcribe())
        with c2: st.write("**Protein Sequence:**"); st.code(s.translate())

# --- 탭 3: 포맷 변환기 ---
with tab3:
    st.subheader("Format Converter")
    raw_data = st.text_area("Raw Text Input", height=150)
    if st.button("Convert to FASTA"):
        if raw_data: st.code(f">Converted_Result\n{raw_data.strip()}")

# --- 탭 4: 프라이머 설계 ---
with tab4:
    st.subheader("Simple Primer Designer")
    target_dna = st.text_area("Target DNA", key="pri_in").upper().strip()
    p_len = st.slider("Primer Length", 18, 25, 20)
    if st.button("Design Primers"):
        if len(target_dna) > p_len*2:
            f_p = Seq(target_dna[:p_len])
            r_p = Seq(target_dna[-p_len:]).reverse_complement()
            st.success(f"Forward: {f_p} (Tm: {mt.Tm_NN(f_p):.1f}℃) | Reverse: {r_p} (Tm: {mt.Tm_NN(r_p):.1f}℃)")
        else: st.warning("서열이 부족합니다.")

# --- 탭 5: 시각화 (DNA, RNA, Protein 선택형) ---
with tab5:
    st.subheader("📊 Multi-Sequence Visualization")
    c_sel, c_in = st.columns([1, 3])
    with c_sel:
        st_type = st.radio("Type", ["DNA", "RNA", "Protein"])
        w_size = st.slider("Window", 5, 50, 20)
        l_color = st.color_picker("Color", "#0984e3")
    with c_in:
        viz_data = st.text_area(f"Input {st_type} for Graph", height=200).upper().strip()
    
    if viz_data and len(viz_data) >= w_size:
        vals = []
        if st_type in ["DNA", "RNA"]:
            ylabel = "GC Content (%)"
            vals = [(viz_data[i:i+w_size].count('G')+viz_data[i:i+w_size].count('C'))/w_size*100 for i in range(len(viz_data)-w_size+1)]
        else:
            ylabel = "Hydrophobicity"
            h_scale = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
            vals = [sum([h_scale.get(aa,0) for aa in viz_data[i:i+w_size]])/w_size for i in range(len(viz_data)-w_size+1)]
        
        if vals:
            fig, ax = plt.subplots(figsize=(10, 4))
            ax.plot(vals, color=l_color, linewidth=2)
            ax.fill_between(range(len(vals)), vals, color=l_color, alpha=0.1)
            ax.set_ylabel(ylabel); ax.set_xlabel("Position"); ax.set_title(f"{st_type} Analysis")
            st.pyplot(fig)

