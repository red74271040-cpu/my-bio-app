import streamlit as st
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
import io
import ssl
import os

# 보안 및 환경 설정
ssl._create_default_https_context = ssl._create_unverified_context
os.environ['CURL_CA_BUNDLE'] = ''

# --- 페이지 설정 ---
st.set_page_config(page_title="Bioinformatics Toolset", layout="wide")

# --- 커스텀 CSS (화이트 & 그레이 미니멀) ---
st.markdown("""
    <style>
    .stApp { background-color: #ffffff; }
    .main-title { font-family: 'Inter', sans-serif; font-size: 28px; font-weight: 700; color: #2D3436; margin-bottom: 20px; }
    .result-card { background-color: #FDFDFD; border: 1px solid #EAEAEA; padding: 20px; border-radius: 4px; margin-bottom: 15px; }
    .stTabs [data-baseweb="tab-list"] { gap: 24px; }
    .stTabs [data-baseweb="tab"] { height: 50px; white-space: pre-wrap; background-color: #F1F2F6; border-radius: 4px 4px 0 0; padding: 10px 20px; }
    .stTabs [aria-selected="true"] { background-color: #2D3436 !important; color: white !important; }
    .guide-box { background-color: #F8F9FA; border-left: 4px solid #D1D4D9; padding: 15px; color: #636E72; font-size: 14px; margin-top: 10px; border-radius: 2px; }
    .target-label { background-color: #F1F2F6; color: #2D3436; padding: 4px 12px; border-radius: 2px; font-weight: 600; font-size: 13px; display: inline-block; margin-bottom: 10px; }
    </style>
    """, unsafe_allow_html=True)

# --- 1. 로그인 체크 로직 ---
if 'authenticated' not in st.session_state:
    st.session_state['authenticated'] = False

def check_password():
    if st.session_state["password_input"] == "knu2026": # 암호 설정
        st.session_state["authenticated"] = True
        del st.session_state["password_input"]
    else:
        st.error("비밀번호가 일치하지 않습니다.")

if not st.session_state['authenticated']:
    st.markdown('<div style="text-align:center; margin-top:150px;">', unsafe_allow_html=True)
    st.title("Restricted Access")
    st.write("본 시스템은 허가된 사용자만 이용 가능합니다.")
    st.text_input("Access Password", type="password", on_change=check_password, key="password_input")
    st.markdown('</div>', unsafe_allow_html=True)
    st.stop()

# --- 2. 기능 함수 정의 ---
def get_pwn_target_analysis(title):
    title = title.lower()
    if any(word in title for word in ["cellulase", "pectate lyase", "eng-1", "pel1", "expansin"]):
        return "식물 세포벽 분해 효소: 소나무 조직 침입 및 파괴 관여"
    elif any(word in title for word in ["acetylcholinesterase", "flp-", "nlp-", "unc-", "motor"]):
        return "신경 및 운동 기관: 선충의 이동성 및 물리적 확산 관여"
    elif any(word in title for word in ["vitellogenin", "cpi-1", "cystatin", "reproduction"]):
        return "생식 및 발달 인자: 선충의 개체수 증식 및 번식 관여"
    elif any(word in title for word in ["v-atpase", "hsp90", "ribosomal", "atp synthase"]):
        return "핵심 생존 유전자: 생명 유지 및 에너지 대사 관여 (고효율 타겟)"
    elif any(word in title for word in ["cytochrome p450", "cyp-", "gst-", "abc transporter"]):
        return "약제 저항성 및 해독: 살선충제 감수성 및 방어 기제 관여"
    elif any(word in title for word in ["tlp", "thaumatin", "vap1", "venom"]):
        return "병원성 분비 단백질: 식물 면역 체계 교란 및 위조 증상 유발"
    else:
        return "기타 유전자 분석: 소나무재선충 특이 서열로 추가 기능 확인 필요"

# --- 3. 메인 레이아웃 ---
col_t, col_l = st.columns([9, 1])
with col_t:
    st.markdown('<p class="main-title">Bioinformatics Analysis Toolset</p>', unsafe_allow_html=True)
with col_l:
    if st.button("Logout"):
        st.session_state['authenticated'] = False
        st.rerun()

tab1, tab2, tab3, tab4 = st.tabs(["🌲 PWN Analysis", "🧬 Central Dogma", "📂 Converter", "🔍 Primer Design"])

# --- 탭 1: 소나무재선충 분석 ---
with tab1:
    st.subheader("B. xylophilus Target Analysis")
    c1_in, c1_gui = st.columns([3, 2])
    with c1_in:
        sequence = st.text_area("DNA Sequence Input", height=150, placeholder="분석할 서열을 입력하십시오.", key="pwn_single")
        if st.button("RUN ANALYSIS", use_container_width=True):
            if not sequence or len(sequence) < 15:
                st.warning("15bp 이상의 서열을 입력해 주십시오.")
            else:
                with st.spinner("NCBI 데이터베이스 검색 중..."):
                    try:
                        res = NCBIWWW.qblast("blastn", "nt", sequence, expect=10, short_query=True, entrez_query="Bursaphelenchus xylophilus [ORGN]")
                        rec = NCBIXML.read(res)
                        if rec.alignments:
                            for i, aln in enumerate(rec.alignments[:5]):
                                st.markdown(f"""<div class="result-card">
                                    <div class="target-label">Candidate {i+1}</div>
                                    <div style="font-weight:700;">{aln.accession}</div>
                                    <div style="font-size:13px; color:#636E72;">{aln.title[:150]}...</div>
                                    <div style="font-size:14px; margin-top:8px;"><b>분석 결과:</b> {get_pwn_target_analysis(aln.title)}</div>
                                </div>""", unsafe_allow_html=True)
                        else: st.info("검출된 상동 서열이 없습니다.")
                    except Exception as e: st.error(f"Error: {e}")
    with c1_gui:
        st.markdown('<div class="guide-box"><b>Analysis Guide</b><br>- 소나무재선충 전용 BLAST 분석 도구입니다.<br>- 상동 유전자의 기능을 카테고리별로 자동 분류합니다.</div>', unsafe_allow_html=True)

# --- 탭 2: 전사 및 번역 ---
with tab2:
    st.subheader("Transcription & Translation")
    dna_in = st.text_area("Enter DNA Sequence", height=100, placeholder="ATGC...", key="dogma_in").strip().upper()
    if dna_in:
        try:
            s = Seq(dna_in)
            cl1, cl2, cl3 = st.columns(3)
            with cl1: st.markdown("**RNA**"); st.code(s.transcribe())
            with cl2: st.markdown("**Protein**"); st.code(s.translate())
            with cl3: st.markdown("**Complementary DNA**"); st.code(s.complement())
        except: st.error("유효하지 않은 서열입니다.")
    st.markdown('<div class="guide-box"><b>Guide:</b> DNA 서열을 RNA 및 단백질 서열로 변환합니다.</div>', unsafe_allow_html=True)

# --- 탭 3: 형식 변환기 ---
with tab3:
    st.subheader("Sequence Format Converter")
    cf1, cf2 = st.columns(2)
    with cf1:
        fi = st.selectbox("From", ["fasta", "genbank", "fastq"])
        rv = st.text_area("Raw Data", height=200, key="conv_in")
    with cf2:
        fo = st.selectbox("To", ["genbank", "fasta", "text"])
        if st.button("CONVERT"):
            if rv:
                try:
                    hi, ho = io.StringIO(rv), io.StringIO()
                    SeqIO.write(list(SeqIO.parse(hi, fi)), ho, fo)
                    st.text_area("Result", value=ho.getvalue(), height=200)
                except Exception as e: st.error(f"Error: {e}")
    st.markdown('<div class="guide-box"><b>Guide:</b> 생물정보학 표준 포맷을 상호 변환합니다.</div>', unsafe_allow_html=True)

# --- 탭 4: 프라이머 설계 ---
with tab4:
    st.subheader("Simple Primer Designer")
    cp1, cp2 = st.columns([3, 2])
    with cp1:
        ts = st.text_area("Enter Target Gene Sequence", height=150, placeholder="ATGC...", key="pri_in").strip().upper()
        p_len = st.slider("Primer Length", 18, 30, 20)
        if st.button("GENERATE PRIMER PAIR"):
            if len(ts) < p_len * 2: st.warning("서열이 너무 짧습니다.")
            else:
                try:
                    f_s, r_s = Seq(ts[:p_len]), Seq(ts[-p_len:]).reverse_complement()
                    for label, s in [("Forward", f_s), ("Reverse", r_s)]:
                        tm = mt.Tm_NN(s)
                        gc = (s.count("G") + s.count("C")) / len(s) * 100
                        st.markdown(f"""<div class="result-card">
                            <div class="target-label">{label} Primer</div>
                            <div style="font-family:monospace; font-size:16px;">{s}</div>
                            <p style="font-size:13px; color:#636E72;">Tm: {tm:.2f}°C | GC: {gc:.1f}%</p>
                        </div>""", unsafe_allow_html=True)
                except Exception as e: st.error(f"Error: {e}")
    with cp2:
        st.markdown('<div class="guide-box"><b>Primer Guide</b><br>- 유전자 양 끝단을 기준으로 프라이머 쌍을 설계합니다.<br>- Reverse는 상보적 역서열로 자동 변환됩니다.</div>', unsafe_allow_html=True)