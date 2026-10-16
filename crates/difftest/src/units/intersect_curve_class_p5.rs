//! Differential runner for the part-5 `SmCurveClassification` operations.
//! Real line/circle sources and line/plane targets; complete interval and solution
//! state is compared. Topology-handle and part-6 gaps are explicit deferrals.

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use solidcore::curve::curve::p2::Curve;
use solidcore::curve::ellipse::{Ellipse, EllipseOptions};
use solidcore::curve::line::{BoundedLineOptions, Line};
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::solution_array::{Solution, SolutionType};
use solidcore::intersect::curve_class::{CurveClassification, CurveInterval, ResetOptions};
use solidcore::intersect::curve_class::p5::{insert_solutions, ExtractOptions};
use solidcore::intersect::point_class::{ClassGeometry,ClassObject,ClassObjectKind,PointClassType,PointClassification};
use solidcore::kernel::Kernel;
use solidcore::math::Vec3;
use solidcore::surface::plane::Plane;
use solidcore::surface::surface::Surface;
use solidcore::surface::volume::Volume;
use solidcore::{Result as ScResult, SolidError};
use solidcore_ffi::units::intersect_curve_class_p5 as ffi;
use crate::case::CaseContext;
use crate::runner::default_cases;
use crate::{Case,HarnessError,Output,Rng,UnitDiff};

pub struct Runner;

fn scale() -> f64 {1.0}
fn zone() -> f64 {1e-6}
fn none_tol() -> f64 {-1.0}
fn yes() -> bool {true}
fn bounds() -> Vec<f64> {vec![0.,1.]}

#[derive(Debug,Clone,Serialize,Deserialize)]
struct Input {
    op: usize,
    #[serde(default)] shape: usize,
    #[serde(default)] dim2: bool,
    #[serde(default)] ignore: bool,
    #[serde(default)] flag1: bool,
    #[serde(default)] flag2: bool,
    #[serde(default="yes")] flag3: bool,
    #[serde(default)] target: usize,
    #[serde(default)] index: usize,
    #[serde(default="scale")] scale: f64,
    #[serde(default="zone")] tol: f64,
    #[serde(default)] angle: f64,
    #[serde(default="bounds")] bounds: Vec<f64>,
    #[serde(default)] other_bounds: Option<Vec<f64>>,
    #[serde(default)] classes: Vec<i32>,
    #[serde(default)] other_classes: Vec<i32>,
    #[serde(default)] solutions: Vec<[f64;7]>,
    #[serde(default="none_tol")] xsect_tol: f64,
}

impl Input {
    fn flat(&self) -> Vec<f64> {
        let n=self.bounds.len().saturating_sub(1);
        let mut out=vec![self.op as f64,self.shape as f64,if self.dim2 {0.} else {1.},f64::from(self.ignore),f64::from(self.flag1),f64::from(self.flag2),f64::from(self.flag3),self.target as f64,self.index as f64,self.scale,self.tol,self.angle,n as f64,self.solutions.len() as f64,self.xsect_tol];
        out.extend(&self.bounds);
        out.extend(self.other_bounds.as_ref().unwrap_or(&self.bounds));
        for c in [&self.classes,&self.other_classes] {
            if c.is_empty() {out.extend(std::iter::repeat_n(100.,3*n));}
            else {out.extend(c.iter().map(|x|f64::from(*x)));}
        }
        out.extend(self.solutions.iter().flatten());out
    }
}

#[derive(Debug,Clone,Copy,PartialEq,Eq)]
struct H(usize,ClassObjectKind);
impl ClassObject for H {fn class_kind(&self)->ClassObjectKind {self.1}}
struct Pool {source:Curve,target:Curve,plane:Surface,vertex:Vec3,vertex_tol:f64}
impl ClassGeometry<H> for Pool {
    fn curve(&self,o:H)->Option<&Curve> {if o.1!=ClassObjectKind::Curve {None} else if o.0==0 {Some(&self.source)} else {Some(&self.target)}}
    fn surface(&self,o:H)->Option<&Surface> {if o.1==ClassObjectKind::Surface {Some(&self.plane)} else {None}}
    fn volume(&self,_o:H)->Option<&Volume> {None}
    fn point(&self,_o:H)->Option<Vec3> {None}
    fn vertex_point(&self,o:H)->Option<Vec3> {if o.1==ClassObjectKind::Vertex {Some(self.vertex)} else {None}}
    fn topology_zone_tol_3d(&self,o:H)->Option<f64> {if o.1==ClassObjectKind::Vertex {Some(self.vertex_tol)} else {None}}
}
fn target(i:&Input)->H {H(1,match i.target {0=>ClassObjectKind::Curve,1=>ClassObjectKind::Surface,2=>ClassObjectKind::Vertex,_=>ClassObjectKind::Face})}

fn setup(k:&Kernel,pool:&Pool,i:&Input,other:bool)->ScResult<CurveClassification<H>> {
    let b=if other {i.other_bounds.as_ref().unwrap_or(&i.bounds)} else {&i.bounds};
    let classes=if other {&i.other_classes} else {&i.classes};
    let mut cc=CurveClassification::new(k);
    cc.reset(k,&ResetOptions {curve:Some(H(0,ClassObjectKind::Curve)),interval:Some(Extent1d::from_min_max(b[0],*b.last().unwrap())),zone_tol_3d:i.tol,..ResetOptions::default()})?;
    cc.set_ignore_tangency(i.ignore);
    cc.set_size(k,b.len()-1)?;
    for (j,bounds) in b.windows(2).enumerate() {
        let mut ivl=CurveInterval::with_interval(k,&Extent1d::from_min_max(bounds[0],bounds[1]),Some(&cc));
        for (w,pc) in [&mut ivl.start,&mut ivl.mid,&mut ivl.end].into_iter().enumerate() {
            let class=PointClassType::from_code(classes.get(3*j+w).copied().unwrap_or(100)).unwrap();
            if class!=PointClassType::Unknown {pc.set_class_object(class,if class==PointClassType::Nothing {None} else {Some(target(i))});}
        }
        *cc.at_mut(j)?=ivl;
    }
    let _=pool; Ok(cc)
}

fn put_pc(out:&mut Vec<f64>,pc:&PointClassification<H>) {
    let uv=pc.uv_param();let uvw=pc.uvw_param();let sss=pc.sss_point();
    out.extend([f64::from(pc.point_class().code()),pc.object().map_or(-1.,|h|if h.0==1 {1.} else {-2.}),pc.deviation(),pc.t_param(),uv.x,uv.y,uvw.x,uvw.y,uvw.z,f64::from(pc.are_parameters_set()),f64::from(pc.sss_point_flag()),sss.x,sss.y,sss.z,pc.pre_snap_param(),pc.src_zone_tol_3d()]);
}
fn put_cc(out:&mut Vec<f64>,cc:&CurveClassification<H>) {
    out.extend([cc.size() as f64,cc.interval().min(),cc.interval().max()]);
    for ivl in cc.intervals() {out.extend([ivl.min(),ivl.max()]);for pc in [&ivl.start,&ivl.mid,&ivl.end] {put_pc(out,pc);}}
}
fn status(e:&SolidError)->f64 {match e {SolidError::Upstream {code,..}=>f64::from(*code),SolidError::InvalidInput {..}=>1007.,_=>-1.}}

fn rust(i:&Input)->ScResult<Output> {
    let k=Kernel::default();
    let source=if i.shape==0 {
        Curve::from(Line::bounded(Vec3::ZERO,Vec3::new(1.,0.,0.),&Extent1d::from_min_max(0.,1.),&BoundedLineOptions {scale:i.scale,dim:if i.dim2 {2} else {3},..BoundedLineOptions::default()}))
    } else {
        Curve::from(Ellipse::new_circle(&k, Vec3::ZERO,Vec3::new(1.,0.,0.),Vec3::new(0.,1.,0.),Extent1d::from_min_max(0.,360.),i.scale,&EllipseOptions {dim:if i.dim2 {2} else {3},..EllipseOptions::default()})?)
    };
    let pool=Pool {source,target:Curve::from(Line::bounded(Vec3::ZERO,Vec3::new(i.angle.cos(),i.angle.sin(),0.),&Extent1d::from_min_max(0.,1.),&BoundedLineOptions {scale:i.scale,dim:if i.dim2 {2} else {3},..BoundedLineOptions::default()})),plane:Surface::from(Plane::from_normal(Vec3::ZERO,Vec3::new(0.,0.,1.))),vertex:Vec3::new(i.angle*i.scale,0.,0.),vertex_tol:i.tol};
    let mut cc=setup(&k,&pool,i,false)?;let mut other=setup(&k,&pool,i,true)?;
    let mut solutions=Vec::new();
    for s in &i.solutions {
        let mut sol=Solution::new();sol.solution_type=SolutionType::from_code(s[0] as i32).unwrap();sol.num_objects=2;sol.num_variables=2;sol.objects[0]=Some(H(0,ClassObjectKind::Curve));sol.objects[1]=Some(target(i));
        *sol.start.param_mut(0)?=s[1];*sol.start.param_mut(1)?=s[2];*sol.end.param_mut(0)?=s[3];*sol.end.param_mut(1)?=s[4];sol.start.solution_value=s[5];sol.end.solution_value=s[6];solutions.push(sol);
    }
    let mut out=Vec::new();
    let result=match i.op {
        0=>insert_solutions(&k,&pool,&mut cc,&mut solutions,target(i),if i.xsect_tol<0. {None} else {Some(i.xsect_tol)}),
        1=>cc.insert_topology_intersections(&k,&pool,&solutions,i.index),
        2=>cc.fixup_classification(&k,&pool),
        3=>cc.extract_geometry(&k,&pool,&mut other,ExtractOptions {uv_curves1:i.flag1,uv_curves2:i.flag2,points_3d:i.flag3}).map(|g|{out.extend([g.curves_3d.len() as f64,g.uv_curves1.len() as f64,g.uv_curves2.len() as f64,g.points_3d.len() as f64]);}),
        4=>cc.merge_classifications(&k,&pool,&mut other,i.flag1).map(|(a,b)|out.extend([a.faces.len() as f64,a.edges.len() as f64,a.vertices.len() as f64,b.faces.len() as f64,b.edges.len() as f64,b.vertices.len() as f64])),
        5=>cc.merge_intervals(i.flag1).map(|a|out.extend([a.faces.len() as f64,a.edges.len() as f64,a.vertices.len() as f64,0.,0.,0.])),
        6=>cc.intersect_vertices(&k,&pool,if i.flag1 {&[H(1,ClassObjectKind::Vertex)]} else {&[]},&pool.plane,i.flag2,i.flag3),
        _=>unreachable!(),
    };
    let st=match result {Ok(())=>0.,Err(e)=>{if let Some(deferred)=Output::unported(&e) {return Ok(deferred);}status(&e)}};
    out.insert(0,st);put_cc(&mut out,&cc);put_cc(&mut out,&other);
    for s in solutions {out.extend([f64::from(s.solution_type.code()),s.start.param(0)?,s.start.param(1)?,s.end.param(0)?,s.end.param(1)?,s.start.solution_value,s.end.solution_value]);}
    Ok(Output::new().values("out",out))
}

impl UnitDiff for Runner {
    fn id(&self)->&'static str {"intersect.curve_class_p5"}
    fn generate(&self,rng:&mut Rng)->Option<Value> {
        for _ in 0..64 {
            let op=rng.below(7) as usize;
            let n=1+rng.below(5) as usize;
            let target=if op==6 {2} else {rng.below(2) as usize};
            let scale=[1e-6,1.,1e6][rng.below(3) as usize];
            let mut i:Input=serde_json::from_value(json!({"op":op,"target":target,"scale":scale})).ok()?;
            i.shape=rng.below(2) as usize;i.dim2=rng.below(4)==0;i.ignore=rng.below(2)==0;i.angle=if rng.below(2)==0 {0.} else {rng.unit()*std::f64::consts::PI};
            i.bounds=(0..=n).map(|j|j as f64/n as f64).collect();
            i.tol=[0.,1e-9,1e-6,.01][rng.below(4) as usize];
            i.flag1=rng.below(2)==0;i.flag2=rng.below(2)==0;i.flag3=rng.below(2)==0;i.index=rng.below(2) as usize;
            if op==0 || op==1 {
                i.xsect_tol=if rng.below(2)==0 {-1.} else {rng.unit()*.01};
                for _ in 0..1+rng.below(4) {
                    let a=rng.unit();let b=if rng.below(2)==0 {a+1e-10} else {rng.unit()};
                    let typ=if rng.below(3)==0 {2.} else {1.};
                    i.solutions.push([typ,a,rng.unit(),b.min(1.),rng.unit(),rng.unit()*1e-7,rng.unit()*1e-7]);
                }
            }
            if ffi::run(&i.flat()).is_ok() {return serde_json::to_value(i).ok();}
        }
        None
    }
    fn cases(&self,ctx:&CaseContext)->Result<Vec<Case>,HarnessError> {
        let context=CaseContext {seeds:ctx.seeds.max(800),..ctx.clone()};
        default_cases(self,&context)
    }
    fn run_cpp(&self,case:&Case)->Result<Output,HarnessError> {let i:Input=case.parse()?;ffi::run(&i.flat()).map(|v|Output::new().values("out",v)).map_err(|e|HarnessError::new(e.to_string()))}
    fn run_rust(&self,case:&Case)->Result<Output,HarnessError> {let i:Input=case.parse()?;rust(&i).map_err(|e|HarnessError::new(e.to_string()))}
}
