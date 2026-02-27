import streamlit as st
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio import SeqIO, Entrez
from Bio.SeqUtils import MeltingTemp as mt
import matplotlib.pyplot as plt
import io, ssl, os

# 보안 및 환경 설정
ssl._create_default_https_context = ssl._create_unverified_context
os.environ['CURL_CA_BUNDLE'] = ''
Entrez.email = "your_email@example.com" # NCBI 에티켓

# --- 페이지 설정 ---
st.set_page_config(page_title="Bio-Research Station", layout="wide")

# --- 커스텀 CSS ---
st.markdown("""
    <style>
    .stApp { background-color: #ffffff; }
    .main-title { font-family: 'Inter', sans-serif; font-size: 28px; font-weight: 700; color: #2D3436; margin-bottom: 20px; }
    .result-card { background-color: #FDFDFD; border: 1px solid #EAEAEA; padding: 20px; border-radius: 4px; margin-bottom: 15px; }
    .target-label { background-color: #F1F2F6; color: #2D3436; padding: 4px 12px; border-radius: 2px; font-weight: 600; font-size: 13px; display: inline-block; margin-bottom: 10px; }
    .guide-box { background-color: #F8F9FA; border-left: 4px solid #D1D4D9; padding: 15px; color: #636E72; font-size: 14px; margin-top: 10px; }
    </style>
    """, unsafe_allow_html=True)

# --- 1. 로그인 체크 로직 ---
if 'auth' not in st.session_state:
    st.session_state['auth'] = False

def check_password():
    if st.session_state["pw_input"] == "knu2026":
        st.session_state['auth'] = True
        del st.session_state["pw_input"]
    else:
        st.error("비밀번호가 일치하지 않습니다.")

if not st.session_state['auth']:
    st.markdown('<div style="text-align:center; margin-top:150px;"><h1>Restricted Access</h1><p>허가된 사용자만 이용 가능합니다.</p></div>', unsafe_allow_html=True)
    st.text_input("Access Password", type="password", on_change=check_password, key="pw_input")
    st.stop()

# --- 2. 분석 함수 정의 ---
def get_pwn_analysis(title):
    t = title.lower()
    if any(w in t for w in ["cellulase", "pectate", "eng-1", "pel1", "expansin"]): return "식물 세포벽 분해 효소: 조직 침입 및 파괴 관여"
    elif any(w in t for w in ["v-atpase", "hsp90", "ribosomal", "atp synthase"]): return "핵심 생존 유전자: 생명 유지 관여 (RNAi 고효율 타겟)"
    elif any(w in t for w in ["cytochrome", "cyp-", "gst-", "abc transporter"]): return "약제 저항성 및 해독: 살선충제 방어 기제 관여"
    return "기타 기능 유전자: 상세 정보를 확인하십시오."

# --- 3. 메인 화면 레이아웃 ---
col_t, col_l = st.columns([9, 1])
with col_t:
    st.markdown('<p class="main-title">Bio-Research Station</p>', unsafe_allow_html=True)
with col_l:
    if st.button("Logout"):
        st.session_state['auth'] = False
        st.rerun()

tab1, tab2, tab3, tab4 = st.tabs(["🌲 Target Analysis", "🧬 Central Dogma", "📂 Converter", "🔍 Primer Design"])

# --- 탭 1: 소나무재선충 분석 & 논문 검색 ---
with tab1:
    st.subheader("B. xylophilus Target & PubMed Search")
    c1_in, c1_gui = st.columns([3, 2])
    with c1_in:
        sequence = st.text_area("DNA Sequence Input", height=150, placeholder="서열을 입력하세요.")
        if st.button("RUN ANALYSIS", use_container_width=True):
            if len(sequence) < 15: st.warning("15bp 이상의 서열이 필요합니다.")
            else:
                with st.spinner("NCBI BLAST & PubMed 검색 중..."):
                    try:
                        res = NCBIWWW.qblast("blastn", "nt", sequence, expect=10, entrez_query="Bursaphelenchus xylophilus [ORGN]")
                        rec = NCBIXML.read(res)
                        if rec.alignments:
                            for i, aln in enumerate(rec.alignments[:3]):
                                st.markdown(f"""<div class="result-card">
                                    <div class="target-label">Candidate {i+1}</div>
                                    <div style="font-weight:700; font-size:18px;">{aln.accession}</div>
                                    <div style="font-size:13px; color:#636E72; margin-bottom:10px;">{aln.title[:150]}...</div>
                                    <div style="font-size:14px;"><b>분석 결과:</b> {get_pwn_analysis(aln.title)}</div>
                                </div>""", unsafe_allow_html=True)
                                
                                # 논문 검색
                                st.write(f"📚 Related Publications for {aln.accession}")
                                try:
                                    sh = Entrez.esearch(db="pubmed", term=f"{aln.accession} Bursaphelenchus", retmax=2)
                                    pids = Entrez.read(sh)["IdList"]
                                    if pids:
                                        for pid in pids:
                                            smh = Entrez.esummary(db="pubmed", id=pid); sm = Entrez.read(smh)
                                            st.markdown(f"- [{sm[0]['Title']}](https://pubmed.ncbi.nlm.nih.gov/{pid}/)")
                                    else: st.write("관련 논문 없음")
                                except: st.write("논문 검색 실패")
                                st.markdown("---")
                        else: st.info("결과가 없습니다.")
                    except Exception as e: st.error(f"Error: {e}")
    with c1_gui:
        st.markdown('<div class="guide-box"><b>Guide</b><br>서열 분석과 동시에 해당 유전자와 관련된 최신 논문 링크를 PubMed에서 가져옵니다.</div>', unsafe_allow_html=True)

# --- 탭 2: 전사/번역 및 GC 그래프 ---
with tab2:
    st.subheader("Transcription & GC Plot")
    d_in = st.text_area("Input DNA Sequence", key="dna_plot").upper().strip()
    if d_in:
        try:
            s = Seq(d_in); cl1, cl2 = st.columns(2)
            with cl1: st.write("**RNA:**"); st.code(s.transcribe())
            with cl2: st.write("**Protein:**"); st.code(s.translate())
            
            # GC 그래프 생성
            win = 20; gcs = [(d_in[i:i+win].count('G')+d_in[i:i+win].count('C'))/win*100 for i in range(len(d_in)-win+1)]
            if gcs:
                fig, ax = plt.subplots(figsize=(10, 3))
                ax.plot(gcs, color='#0984e3'); ax.set_ylim(0, 100); ax.set_ylabel("GC %")
                st.pyplot(fig)
        except: st.error("유효하지 않은 서열입니다.")

# --- 탭 3: 포맷 변환기 ---
with tab3:
    st.subheader("Format Converter")
    cv_in = st.text_area("Raw Data", height=150)
    if st.button("Convert to FASTA"):
        st.code(">Converted_Sequence\n" + cv_in.strip())

# --- 탭 4: 프라이머 설계 ---
with tab4:
    st.subheader("Primer Designer")
    ts = st.text_area("Target Gene", key="pri_target").strip().upper()
    pl = st.slider("Length", 18, 25, 20)
    if st.button("Generate"):
        if len(ts) > pl*2:
            f, r = Seq(ts[:pl]), Seq(ts[-pl:]).reverse_complement()
            for lab, p in [("Forward", f), ("Reverse", r)]:
                tm = mt.Tm_NN(p); gc = (p.count('G')+p.count('C'))/len(p)*100
                st.markdown(f"""<div class="result-card"><b>{lab}:</b> {p}<br>Tm: {tm:.1f}°C | GC: {gc:.1f}%</div>""", unsafe_allow_html=True)
