/*
 * C# / XNA  port of Bullet (c) 2011 Mark Neale <xexuxjy@hotmail.com>
 *
 * Bullet Continuous Collision Detection and Physics Library
 * Copyright (c) 2003-2008 Erwin Coumans  http://www.bulletphysics.com/
 *
 * This software is provided 'as-is', without any express or implied warranty.
 * In no event will the authors be held liable for any damages arising from
 * the use of this software.
 * 
 * Permission is granted to anyone to use this software for any purpose, 
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 * 
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 */

using System;
using System.Collections.Generic;
using System.Diagnostics;
using BulletXNA.BulletCollision;
using BulletXNA.LinearMath;

namespace BulletXNA.BulletDynamics
{
    public class RigidBody : CollisionObject
    {
        private const float _maxAngvel = MathUtil.SIMD_HALF_PI;
        public static int UniqueId = 0;

        private IndexedBasisMatrix	_mInvInertiaTensorWorld= IndexedBasisMatrix.Identity;
	    private IndexedVector3		_mLinearVelocity;
	    private IndexedVector3		_mAngularVelocity;
	    private float		_mInverseMass;
        private IndexedVector3     _mLinearFactor;


	    private IndexedVector3		_mGravity;	
	    private IndexedVector3		_mGravityAcceleration;
	    private IndexedVector3		_mInvInertiaLocal;
	    private IndexedVector3		_mTotalForce;
	    private IndexedVector3		_mTotalTorque;
    	
	    private float		_mLinearDamping;
	    private float		_mAngularDamping;

	    private bool		_mAdditionalDamping;
	    private float		_mAdditionalDampingFactor;
	    private float		_mAdditionalLinearDampingThresholdSqr;
	    private float		_mAdditionalAngularDampingThresholdSqr;
	    private float		_mAdditionalAngularDampingFactor;

	    private float		_mLinearSleepingThreshold;
	    private float		_mAngularSleepingThreshold;

	    //m_optionalMotionState allows to automatic synchronize the world transform for active objects
	    private IMotionState	_mOptionalMotionState;

	    //keep track of typed constraints referencing this rigid body
	    private IList<TypedConstraint> _mConstraintRefs;

        private RigidBodyFlags _mRigidbodyFlags;

        public int MDebugBodyId;

        public IndexedVector3 MDeltaLinearVelocity;
        public IndexedVector3 MDeltaAngularVelocity;
        protected IndexedVector3 MAngularFactor;
        public IndexedVector3 MInvMass;
        protected IndexedVector3 MPushVelocity;
        protected IndexedVector3 MTurnVelocity;

		//static RigidBody()
		//{
		//    String filename = @"C:\users\man\xna-rb-output.txt";
		//    s_filestream = File.Open(filename, FileMode.Create, FileAccess.Write, FileShare.None);
		//    s_streamWriter = new StreamWriter(s_filestream);
		//}

        public RigidBody()
        { }

	    ///btRigidBody constructor using construction info
	    public RigidBody(RigidBodyConstructionInfo constructionInfo)
        {
            SetupRigidBody(constructionInfo);
        }

	    ///btRigidBody constructor for backwards compatibility. 
	    ///To specify friction (etc) during rigid body construction, please use the other constructor (using btRigidBodyConstructionInfo)
	    public RigidBody(float mass, IMotionState motionState, CollisionShape collisionShape, IndexedVector3 localInertia)
        {
            RigidBodyConstructionInfo cinfo = new RigidBodyConstructionInfo(mass,motionState,collisionShape,localInertia);
	        SetupRigidBody(cinfo);
        }

        public override void Cleanup()
        {
            base.Cleanup();
            //No constraints should point to this rigidbody
            //Remove constraints from the dynamics world before you delete the related rigidbodies. 
            Debug.Assert(_mConstraintRefs.Count == 0);

        }

	    ///setupRigidBody is only used internally by the constructor
	    protected void	SetupRigidBody(RigidBodyConstructionInfo constructionInfo)
        {
	        m_internalType=CollisionObjectTypes.CO_RIGID_BODY;

	        _mLinearVelocity = IndexedVector3.Zero;
	        _mAngularVelocity = IndexedVector3.Zero;
            MAngularFactor = IndexedVector3.One;
            _mLinearFactor = IndexedVector3.One;
	        _mGravity = IndexedVector3.Zero;
	        _mGravityAcceleration = IndexedVector3.Zero;
	        _mTotalForce = IndexedVector3.Zero;
	        _mTotalTorque = IndexedVector3.Zero;
			SetDamping(constructionInfo.m_linearDamping, constructionInfo.m_angularDamping);
	        _mLinearSleepingThreshold = constructionInfo.m_linearSleepingThreshold;
	        _mAngularSleepingThreshold = constructionInfo.m_angularSleepingThreshold;
	        _mOptionalMotionState = constructionInfo.m_motionState;
	        m_contactSolverType = 0;
	        m_frictionSolverType = 0;
	        _mAdditionalDamping = constructionInfo.m_additionalDamping;
	        _mAdditionalDampingFactor = constructionInfo.m_additionalDampingFactor;
	        _mAdditionalLinearDampingThresholdSqr = constructionInfo.m_additionalLinearDampingThresholdSqr;
	        _mAdditionalAngularDampingThresholdSqr = constructionInfo.m_additionalAngularDampingThresholdSqr;
	        _mAdditionalAngularDampingFactor = constructionInfo.m_additionalAngularDampingFactor;

	        if (_mOptionalMotionState != null)
	        {
		        _mOptionalMotionState.GetWorldTransform(out m_worldTransform);
	        } 
            else
	        {
		        SetWorldTransform(ref constructionInfo.m_startWorldTransform);
	        }

	        m_interpolationWorldTransform = m_worldTransform;
            m_interpolationLinearVelocity = IndexedVector3.Zero;
            m_interpolationAngularVelocity = IndexedVector3.Zero;
        	
	        //moved to btCollisionObject
	        m_friction = constructionInfo.m_friction;
	        m_restitution = constructionInfo.m_restitution;

	        SetCollisionShape( constructionInfo.m_collisionShape );
	        MDebugBodyId = UniqueId++;
        	
	        SetMassProps(constructionInfo.m_mass, constructionInfo.m_localInertia);
	        UpdateInertiaTensor();
            _mRigidbodyFlags = RigidBodyFlags.BT_NONE;
            _mConstraintRefs = new List<TypedConstraint>();

            MDeltaLinearVelocity = IndexedVector3.Zero;
            MDeltaAngularVelocity = IndexedVector3.Zero;
            MInvMass = _mInverseMass * _mLinearFactor;
            MPushVelocity = IndexedVector3.Zero;
            MTurnVelocity = IndexedVector3.Zero;

        }

        public void ProceedToTransform(ref IndexedMatrix newTrans)
        {
            SetCenterOfMassTransform(ref newTrans);
        }
	
	    ///to keep collision detection and dynamics separate we don't store a rigidbody pointer
	    ///but a rigidbody is derived from btCollisionObject, so we can safely perform an upcast
	    public static RigidBody	Upcast(CollisionObject colObj)
	    {
		    if ((colObj.GetInternalType()&CollisionObjectTypes.CO_RIGID_BODY) != 0)
            {
			    return colObj as RigidBody;
            }
		    return null;
	    }

	    /// continuous collision detection needs prediction
	    public void	PredictIntegratedTransform(float timeStep, out IndexedMatrix predictedTransform) 
        {
#if DEBUG        
			if (BulletGlobals.g_streamWriter != null && BulletGlobals.debugRigidBody)
			{
                BulletGlobals.g_streamWriter.WriteLine("[{0}] predictIntegratedTransform pre", (String)m_userObjectPointer);
				MathUtil.PrintMatrix(BulletGlobals.g_streamWriter,m_worldTransform);
				MathUtil.PrintVector3(BulletGlobals.g_streamWriter,"LinearVel", _mLinearVelocity);
				MathUtil.PrintVector3(BulletGlobals.g_streamWriter,"AngularVel",_mAngularVelocity);
			}
#endif			
            TransformUtil.IntegrateTransform(ref m_worldTransform, ref _mLinearVelocity, ref _mAngularVelocity, timeStep, out predictedTransform);
            MathUtil.SanityCheckVector(m_worldTransform._basis[1]);

#if DEBUG            
			if (BulletGlobals.g_streamWriter != null && BulletGlobals.debugRigidBody)
			{
				BulletGlobals.g_streamWriter.WriteLine("[{0}] predictIntegratedTransform post", (String)m_userObjectPointer);
				MathUtil.PrintMatrix(BulletGlobals.g_streamWriter, predictedTransform);
			}
#endif			
        }

        public void SaveKinematicState(float timeStep)
        {
            //todo: clamp to some (user definable) safe minimum timestep, to limit maximum angular/linear velocities
            if (timeStep != 0f)
            {
                //if we use motionstate to synchronize world transforms, get the new kinematic/animated world transform
                if (GetMotionState() != null)
                {
                    GetMotionState().GetWorldTransform(out m_worldTransform);
                }
                IndexedVector3 linVel = IndexedVector3.Zero, angVel = IndexedVector3.Zero;

                // debug steps to track NaN's
                IndexedMatrix worldTransform = m_worldTransform;
                TransformUtil.CalculateVelocity(ref m_interpolationWorldTransform, ref worldTransform, timeStep, out _mLinearVelocity, out _mAngularVelocity);
                SetWorldTransform(ref worldTransform);

                m_interpolationLinearVelocity = _mLinearVelocity;
                m_interpolationAngularVelocity = _mAngularVelocity;
                SetInterpolationWorldTransform(ref m_worldTransform);
                //printf("angular = %f %f %f\n",m_angularVelocity.getX(),m_angularVelocity.getY(),m_angularVelocity.getZ());
            }
        }

        public void ApplyGravity()
        {
            if (IsStaticOrKinematicObject())
            {
                return;
            }

            ApplyCentralForce(ref _mGravity);
        }

        public void SetGravity(IndexedVector3 acceleration)
        {
            SetGravity(ref acceleration);
        }

        public void SetGravity(ref IndexedVector3 acceleration)
        {
            if (_mInverseMass != 0f)
            {
                _mGravity = acceleration * (1f / _mInverseMass);
            }
            _mGravityAcceleration = acceleration;
        }

        public IndexedVector3 Gravity{
            get { return _mGravityAcceleration; }
        }

        public void SetDamping(float lin_damping, float ang_damping)
        {
            _mLinearDamping = MathUtil.Clamp(lin_damping, 0f, 1f);
            _mAngularDamping = MathUtil.Clamp(ang_damping, 0f, 1f);

        }

	    public float GetLinearDamping()
	    {
		    return _mLinearDamping;
	    }

	    public float GetAngularDamping()
	    {
		    return _mAngularDamping;
	    }

	    public float GetLinearSleepingThreshold()
	    {
		    return _mLinearSleepingThreshold;
	    }

	    public float GetAngularSleepingThreshold() 
	    {
		    return _mAngularSleepingThreshold;
	    }

	    public void	ApplyDamping(float timeStep)
        {
	        //On new damping: see discussion/issue report here: http://code.google.com/p/bullet/issues/detail?id=74
	        //todo: do some performance comparisons (but other parts of the engine are probably bottleneck anyway

        //#define USE_OLD_DAMPING_METHOD 1
        #if USE_OLD_DAMPING_METHOD
	        m_linearVelocity *= GEN_clamped((float(1.) - timeStep * m_linearDamping), (float)float(0.0), (float)float(1.0));
	        m_angularVelocity *= GEN_clamped((float(1.) - timeStep * m_angularDamping), (float)float(0.0), (float)float(1.0));
        #else
	        _mLinearVelocity *= (float)Math.Pow((1f-_mLinearDamping), timeStep);
            _mAngularVelocity *= (float)Math.Pow((1f - _mAngularDamping), timeStep);
            MathUtil.SanityCheckVector(ref _mLinearVelocity);
            MathUtil.SanityCheckVector(ref _mAngularVelocity);
#endif

	        if (_mAdditionalDamping)
	        {
		        //Additional damping can help avoiding lowpass jitter motion, help stability for ragdolls etc.
		        //Such damping is undesirable, so once the overall simulation quality of the rigid body dynamics system has improved, this should become obsolete
		        if ((_mAngularVelocity.LengthSquared() < _mAdditionalAngularDampingThresholdSqr) &&
			        (_mLinearVelocity.LengthSquared() < _mAdditionalLinearDampingThresholdSqr))
		        {
			        _mAngularVelocity *= _mAdditionalDampingFactor;
			        _mLinearVelocity *= _mAdditionalDampingFactor;
		        }


                MathUtil.SanityCheckVector(ref _mLinearVelocity);
                MathUtil.SanityCheckVector(ref _mAngularVelocity);
                
                float speed = _mLinearVelocity.Length();
		        if (speed < _mLinearDamping)
		        {
			        float dampVel = 0.005f;
			        if (speed > dampVel)
			        {
				        IndexedVector3 dir = _mLinearVelocity;
                        dir.Normalize();
				        _mLinearVelocity -=  dir * dampVel;
			        } 
                    else
			        {
				        _mLinearVelocity = IndexedVector3.Zero;
			        }
		        }

		        float angSpeed = _mAngularVelocity.Length();
		        if (angSpeed < _mAngularDamping)
		        {
			        float angDampVel = 0.005f;
			        if (angSpeed > angDampVel)
			        {
				        IndexedVector3 dir = _mAngularVelocity;
                        dir.Normalize();
				        _mAngularVelocity -=  dir * angDampVel;
			        } else
			        {
                        _mAngularVelocity = IndexedVector3.Zero;
			        }
		        }
	        }
            MathUtil.SanityCheckVector(ref _mLinearVelocity);
            MathUtil.SanityCheckVector(ref _mAngularVelocity);

        }

        public void SetMassProps(float mass, IndexedVector3 inertia)
        {
            SetMassProps(mass, ref inertia);
        }

	    public void	SetMassProps(float mass, ref IndexedVector3 inertia)
        {
	        if (MathUtil.FuzzyZero(mass))
	        {
		        m_collisionFlags |= CollisionFlags.CF_STATIC_OBJECT;
		        _mInverseMass = 0f;
	        } 
            else
	        {
		        m_collisionFlags &= (~CollisionFlags.CF_STATIC_OBJECT);
		        _mInverseMass = 1.0f / mass;
	        }

			_mGravity = mass * _mGravityAcceleration;

            _mInvInertiaLocal = new IndexedVector3(
                            (inertia.X != 0f) ? 1f / inertia.X : 0f,
                           (inertia.Y !=  0f) ? 1f / inertia.Y : 0f,
                           (inertia.Z !=  0f) ? 1f / inertia.Z : 0f);
            MInvMass = _mLinearFactor * _mInverseMass;
        }
	
        public IndexedVector3 GetLinearFactor()
	    {
		    return _mLinearFactor;
	    }

        public void SetLinearFactor(IndexedVector3 linearFactor)
	    {
		    _mLinearFactor = linearFactor;
		    MInvMass = _mLinearFactor*_mInverseMass;
	    }

        public float InvMass{
            get { return _mInverseMass; }
        }

        public IndexedBasisMatrix GetInvInertiaTensorWorld()
        { 
		    return _mInvInertiaTensorWorld; 
	    }

        public void IntegrateVelocities(float step)
        {
	        if (IsStaticOrKinematicObject())
		        return;

#if DEBUG
			if (BulletGlobals.g_streamWriter != null && BulletGlobals.debugRigidBody)
			{
                BulletGlobals.g_streamWriter.WriteLine(String.Format("[{0}] RigidBody integrateVelocities", (String)m_userObjectPointer));
				MathUtil.PrintVector3(BulletGlobals.g_streamWriter, "integrate LinVel pre", _mLinearVelocity);
				MathUtil.PrintVector3(BulletGlobals.g_streamWriter, "integrate AngVel pre", _mAngularVelocity);
			}
#endif


	        _mLinearVelocity += _mTotalForce * (_mInverseMass * step);
            MathUtil.SanityCheckVector(ref _mLinearVelocity);
            _mAngularVelocity += _mInvInertiaTensorWorld * _mTotalTorque * step;
            MathUtil.SanityCheckVector(ref _mAngularVelocity);
        
	        /// clamp angular velocity. collision calculations will fail on higher angular velocities	
	        float angvel = _mAngularVelocity.Length();
	        if (angvel*step > _maxAngvel)
	        {
		        _mAngularVelocity *= (_maxAngvel/step) /angvel;
	        }
            MathUtil.SanityCheckVector(ref _mAngularVelocity);

#if DEBUG
			if (BulletGlobals.g_streamWriter != null && BulletGlobals.debugRigidBody)
			{
				MathUtil.PrintVector3(BulletGlobals.g_streamWriter, "integrate LinVel post", _mLinearVelocity);
				MathUtil.PrintVector3(BulletGlobals.g_streamWriter, "integrate AngVel post", _mAngularVelocity);
			}
#endif			
        }

        public void SetCenterOfMassTransform(ref IndexedMatrix xform)
        {
#if DEBUG        
			if (BulletGlobals.g_streamWriter != null && BulletGlobals.debugRigidBody)
			{
				BulletGlobals.g_streamWriter.WriteLine(String.Format("[{0}] RigidBody setCenterOfMassTransform",(String)m_userObjectPointer));
				MathUtil.PrintMatrix(BulletGlobals.g_streamWriter, xform);
			}
#endif

            if (IsKinematicObject())
            {
                SetInterpolationWorldTransform(ref m_worldTransform);
            }
            else
            {
                SetInterpolationWorldTransform(ref xform);
            }
            m_interpolationLinearVelocity = GetLinearVelocity();
            m_interpolationAngularVelocity = GetAngularVelocity();
            SetWorldTransform(ref xform);
            UpdateInertiaTensor();
#if DEBUG
			if (BulletGlobals.g_streamWriter != null && BulletGlobals.debugRigidBody)
			{
				BulletGlobals.g_streamWriter.WriteLine("RigidBody setCenterOfMassTransform after calcs");
				MathUtil.PrintMatrix(BulletGlobals.g_streamWriter, m_worldTransform);
			}
#endif
        }

	    public void ApplyCentralForce(ref IndexedVector3 force)
	    {
            _mTotalForce += force * _mLinearFactor;
        }

	    public IndexedVector3 GetTotalForce()
	    {
		    return _mTotalForce;
	    }

	    public IndexedVector3 GetTotalTorque()
	    {
		    return _mTotalTorque;
	    }
    
	    public IndexedVector3 GetInvInertiaDiagLocal()
	    {
		    return _mInvInertiaLocal;
	    }

	    public void SetInvInertiaDiagLocal(ref IndexedVector3 diagInvInertia)
	    {
		    _mInvInertiaLocal = diagInvInertia;
	    }

	    public void	SetSleepingThresholds(float linear,float angular)
	    {
		    _mLinearSleepingThreshold = linear;
		    _mAngularSleepingThreshold = angular;
	    }

        public void ApplyTorque(IndexedVector3 torque)
        {
            ApplyTorque(ref torque);
        }

	    public void	ApplyTorque(ref IndexedVector3 torque)
	    {
            _mTotalTorque += torque * MAngularFactor;
        }
	
	    public void	ApplyForce(ref IndexedVector3 force, ref IndexedVector3 rel_pos) 
	    {
            ApplyCentralForce(ref force);
            IndexedVector3 tempTorque = IndexedVector3.Cross(rel_pos,force);
            tempTorque *= MAngularFactor;
            ApplyTorque(IndexedVector3.Cross(rel_pos,(force * _mLinearFactor)));
        }
	
	    public void ApplyCentralImpulse(ref IndexedVector3 impulse)
	    {

            _mLinearVelocity += impulse * _mLinearFactor * _mInverseMass;
            MathUtil.SanityCheckVector(ref _mLinearVelocity);
	    }

        public void ApplyTorqueImpulse(IndexedVector3 torque)
        {
            ApplyTorqueImpulse(ref torque);
        }

  	    public void ApplyTorqueImpulse(ref IndexedVector3 torque)
	    {
            _mAngularVelocity += _mInvInertiaTensorWorld * torque * MAngularFactor;
        }

        public void ApplyImpulse(IndexedVector3 impulse, IndexedVector3 rel_pos)
        {
            ApplyImpulse(ref impulse, ref rel_pos);
        }
	
	    public void ApplyImpulse(ref IndexedVector3 impulse, ref IndexedVector3 rel_pos) 
	    {
		    if (_mInverseMass != 0f)
		    {
			    ApplyCentralImpulse(ref impulse);
			    if (MAngularFactor.LengthSquared() > 0f)
			    {
				    ApplyTorqueImpulse(IndexedVector3.Cross(rel_pos,(impulse*_mLinearFactor)));
			    }
		    }
	    }

	    //Optimization for the iterative solver: avoid calculating constant terms involving inertia, normal, relative position
        public void InternalApplyImpulse(IndexedVector3 linearComponent, IndexedVector3 angularComponent, float impulseMagnitude,String caller)
        {
            if (impulseMagnitude > 20f)
            {
                int ibreak = 0;
            }
            InternalApplyImpulse(ref linearComponent, ref angularComponent, impulseMagnitude,caller);
        }
	
	    public void ClearForces() 
	    {
		    _mTotalForce = IndexedVector3.Zero;
		    _mTotalTorque = IndexedVector3.Zero;
	    }
	
	    public void UpdateInertiaTensor()
        {
#if DEBUG        
			if (BulletGlobals.g_streamWriter != null && BulletGlobals.debugRigidBody)
            {
                BulletGlobals.g_streamWriter.WriteLine(String.Format("[{0}] RigidBody updateInertiaTensor",(String)m_userObjectPointer));
                MathUtil.PrintVector3(BulletGlobals.g_streamWriter, "invInertiaLocal", _mInvInertiaLocal);
                MathUtil.PrintMatrix(BulletGlobals.g_streamWriter, m_worldTransform);
                MathUtil.PrintMatrix(BulletGlobals.g_streamWriter, m_worldTransform._basis.Scaled(ref _mInvInertiaLocal));
                MathUtil.PrintMatrix(BulletGlobals.g_streamWriter, m_worldTransform._basis.Transpose());
            }
#endif            
            _mInvInertiaTensorWorld = m_worldTransform._basis.Scaled(ref _mInvInertiaLocal) * m_worldTransform._basis.Transpose();
#if DEBUG
			if (BulletGlobals.g_streamWriter != null && BulletGlobals.debugRigidBody)
            {
                MathUtil.PrintMatrix(BulletGlobals.g_streamWriter,_mInvInertiaTensorWorld);
            }
#endif
        }
	
	    public IndexedVector3 GetCenterOfMassPosition() 
        { 
		    return m_worldTransform._origin; 
	    }
	    
        public IndexedQuaternion GetOrientation()
        {
            return m_worldTransform._basis.GetRotation();
        }
	
	    public IndexedMatrix GetCenterOfMassTransform() 
        { 
		    return m_worldTransform; 
	    }

        public IndexedVector3 GetLinearVelocity()
        {
		    return _mLinearVelocity; 
	    }

	    public IndexedVector3 GetAngularVelocity() 
        { 
		    return _mAngularVelocity; 
	    }

        public void SetLinearVelocity(IndexedVector3 lin_vel)
        {
            SetLinearVelocity(ref lin_vel);
        }

	    public void SetLinearVelocity(ref IndexedVector3 lin_vel)
	    { 
		    _mLinearVelocity = lin_vel;
            MathUtil.SanityCheckVector(ref _mLinearVelocity);
        }

        public void SetAngularVelocity(IndexedVector3 ang_vel)
        {
            SetAngularVelocity(ref ang_vel);
        }

	    public void SetAngularVelocity(ref IndexedVector3 ang_vel) 
	    { 
		    _mAngularVelocity = ang_vel; 
	    }

	    public IndexedVector3 GetVelocityInLocalPoint(ref IndexedVector3 rel_pos)
	    {
		    //we also calculate lin/ang velocity for kinematic objects

            IndexedVector3 temp = new IndexedVector3(_mAngularVelocity.Y * rel_pos.Z - _mAngularVelocity.Z * rel_pos.Y,
                _mAngularVelocity.Z * rel_pos.X - _mAngularVelocity.X * rel_pos.Z,
                _mAngularVelocity.X * rel_pos.Y - _mAngularVelocity.Y * rel_pos.X);

            return new IndexedVector3(_mLinearVelocity.X + temp.X, _mLinearVelocity.Y + temp.Y, _mLinearVelocity.Z + temp.Z);

            //return m_linearVelocity + IndexedVector3.Cross(m_angularVelocity,rel_pos);

		    //for kinematic objects, we could also use use:
		    //		return 	(m_worldTransform(rel_pos) - m_interpolationWorldTransform(rel_pos)) / m_kinematicTimeStep;
	    }

	    public void Translate(ref IndexedVector3 v) 
	    {
		    m_worldTransform._origin += v; 
	    }
	
	    public void	GetAabb(out IndexedVector3 aabbMin,out IndexedVector3 aabbMax)
        {
            GetCollisionShape().GetAabb(m_worldTransform, out aabbMin, out aabbMax);
        }
	
	    public float ComputeImpulseDenominator(ref IndexedVector3 pos, ref IndexedVector3 normal)
	    {
		    IndexedVector3 r0 = pos - GetCenterOfMassPosition();

            IndexedVector3 c0 = r0.Cross(ref normal);

            IndexedVector3 vec = (c0 * GetInvInertiaTensorWorld()).Cross(ref r0);

		    return _mInverseMass + IndexedVector3.Dot(normal,vec);

	    }

	    public float ComputeAngularImpulseDenominator(ref IndexedVector3 axis)
	    {
            IndexedVector3 vec = axis * GetInvInertiaTensorWorld();
            return axis.Dot(ref vec);
        }

	    public void	UpdateDeactivation(float timeStep)
	    {
            if ((GetActivationState() == ActivationState.ISLAND_SLEEPING) || (GetActivationState() == ActivationState.DISABLE_DEACTIVATION))
            {
			    return;
            }
		    
            if ((GetLinearVelocity().LengthSquared() < _mLinearSleepingThreshold*_mLinearSleepingThreshold) &&
			    (GetAngularVelocity().LengthSquared() < _mAngularSleepingThreshold*_mAngularSleepingThreshold))
		    {
			    m_deactivationTime += timeStep;
		    } 
            else
		    {
			    m_deactivationTime=0f;
			    SetActivationState(ActivationState.UNDEFINED);
		    }

	    }

	    public bool	WantsSleeping()
	    {

		    if (GetActivationState() == ActivationState.DISABLE_DEACTIVATION)
            {
			    return false;
            }

		    //disable deactivation
            if (BulletGlobals.gDisableDeactivation || BulletGlobals.gDeactivationTime == 0f)
            {
			    return false;
            }

            if ((GetActivationState() == ActivationState.ISLAND_SLEEPING) || (GetActivationState() == ActivationState.WANTS_DEACTIVATION))
            {
			    return true;
            }

            if (m_deactivationTime > BulletGlobals.gDeactivationTime)
		    {
			    return true;
		    }
		    return false;
	    }
	
	    public BroadphaseProxy	GetBroadphaseProxy() 
	    {
		    return m_broadphaseHandle;
	    }

	    public void	SetNewBroadphaseProxy(BroadphaseProxy broadphaseProxy)
	    {
		    m_broadphaseHandle = broadphaseProxy;
	    }

	    //btMotionState allows to automatic synchronize the world transform for active objects
	    public IMotionState	GetMotionState()
	    {
		    return _mOptionalMotionState;
	    }

	    public void	SetMotionState(IMotionState motionState)
	    {
		    _mOptionalMotionState = motionState;
		    if (_mOptionalMotionState != null)
            {
                motionState.GetWorldTransform(out m_worldTransform);
            }
	    }

	    //for experimental overriding of friction/contact solver func
	    int	m_contactSolverType;
	    int	m_frictionSolverType;

	    public void	SetAngularFactor(float angFac)
	    {
		    MAngularFactor = new IndexedVector3(angFac);
	    }

        public void SetAngularFactor(IndexedVector3 angFac)
		{
			SetAngularFactor(ref angFac);
		}

        public void SetAngularFactor(ref IndexedVector3 angFac)
	    {
		    MAngularFactor = angFac;
	    }

	    public IndexedVector3 GetAngularFactor()
	    {
		    return MAngularFactor;
	    }

        public void	SetFlags(RigidBodyFlags flags)
	    {
		    _mRigidbodyFlags = flags;
	    }

	    public RigidBodyFlags GetFlags()
	    {
		    return _mRigidbodyFlags;
	    }

	public IndexedVector3 GetDeltaLinearVelocity()
	{
		return MDeltaLinearVelocity;
	}

	public IndexedVector3 GetDeltaAngularVelocity() 
	{
		return MDeltaAngularVelocity;
	}

	public IndexedVector3 GetPushVelocity()
	{
		return MPushVelocity;
	}

	public IndexedVector3 GetTurnVelocity() 
	{
		return MTurnVelocity;
	}


	    //is this rigidbody added to a btCollisionWorld/btDynamicsWorld/btBroadphase?
	    public bool IsInWorld()
	    {
		    return (GetBroadphaseProxy() != null);
	    }

        public override bool CheckCollideWithOverride(CollisionObject co)
        {
	        RigidBody otherRb = RigidBody.Upcast(co);
	        if (otherRb == null)
		        return true;

	        for (int i = 0; i < _mConstraintRefs.Count; ++i)
	        {
		        TypedConstraint c = _mConstraintRefs[i];
                if (c.IsEnabled())
                {
                    if (c.GetRigidBodyA() == otherRb || c.GetRigidBodyB() == otherRb)
                    {
                        return false;
                    }
                }
	        }

	        return true;

        }

        public void AddConstraintRef(TypedConstraint c)
        {
            if (!_mConstraintRefs.Contains(c))
            {
                _mConstraintRefs.Add(c);
            }

            m_checkCollideWith = true;
        }
        public void RemoveConstraintRef(TypedConstraint c)
        {
            _mConstraintRefs.Remove(c);
            m_checkCollideWith = _mConstraintRefs.Count > 0;

        }

	    public TypedConstraint GetConstraintRef(int index)
	    {
		    return _mConstraintRefs[index];
	    }

	    public int GetNumConstraintRefs()
	    {
		    return _mConstraintRefs.Count;
	    }

        	////////////////////////////////////////////////
	    ///some internal methods, don't use them
    		
	    public IndexedVector3 InternalGetDeltaLinearVelocity()
	    {
		    return MDeltaLinearVelocity;
	    }

        public void InternalSetDeltaLinearVelocity(ref IndexedVector3 v)
        {
            MDeltaLinearVelocity = v;
            MathUtil.SanityCheckVector(ref MDeltaLinearVelocity);
        }

	    public IndexedVector3 InternalGetDeltaAngularVelocity()
	    {
		    return MDeltaAngularVelocity;
        }

        public void InternalSetDeltaAngularVelocity(ref IndexedVector3 v)
        {
            MDeltaAngularVelocity = v;
            MathUtil.SanityCheckVector(ref MDeltaAngularVelocity);
        }

	    public IndexedVector3 InternalGetAngularFactor()
	    {
		    return MAngularFactor;
	    }

	    public IndexedVector3 InternalInvMass {
            get { return MInvMass; }
	    }
    	
	    public IndexedVector3 InternalGetPushVelocity()
	    {
		    return MPushVelocity;
	    }

	    public IndexedVector3 InternalGetTurnVelocity()
	    {
		    return MTurnVelocity;
	    }

        public void InternalSetTurnVelocity(ref IndexedVector3 velocity)
        {
            MTurnVelocity = velocity;
        }

        public void InternalSetPushVelocity(ref IndexedVector3 velocity)
        {
            MPushVelocity = velocity;
        }


	    public void	InternalGetVelocityInLocalPointObsolete(ref IndexedVector3 rel_pos, ref IndexedVector3 velocity )
	    {
		    velocity = GetLinearVelocity()+MDeltaLinearVelocity + IndexedVector3.Cross((GetAngularVelocity()+MDeltaAngularVelocity),rel_pos);
	    }

	    public void	InternalGetAngularVelocity(ref IndexedVector3 angVel)
	    {
		    angVel = GetAngularVelocity()+MDeltaAngularVelocity;
	    }

	    //Optimization for the iterative solver: avoid calculating constant terms involving inertia, normal, relative position
        public void InternalApplyImpulse(ref IndexedVector3 linearComponent, ref IndexedVector3 angularComponent, float impulseMagnitude, String caller)
	    {
#if DEBUG	    
			if (BulletGlobals.g_streamWriter != null && BulletGlobals.debugRigidBody)
			{
				BulletGlobals.g_streamWriter.WriteLine(String.Format("[{0}] internalApplyImpule [{1}]", (String)m_userObjectPointer,caller));
				MathUtil.PrintVector3(BulletGlobals.g_streamWriter, "linComponenet", linearComponent);
				MathUtil.PrintVector3(BulletGlobals.g_streamWriter, "angComponenet", angularComponent);
				BulletGlobals.g_streamWriter.WriteLine("magnitude [{0:0.00000000}]", impulseMagnitude);
			}
#endif
		    if (_mInverseMass != 0f)
            {   
                MDeltaLinearVelocity.X += impulseMagnitude * linearComponent.X;
                MDeltaLinearVelocity.Y += impulseMagnitude * linearComponent.Y;
                MDeltaLinearVelocity.Z += impulseMagnitude * linearComponent.Z;
                //m_deltaLinearVelocity += linearComponent*impulseMagnitude;
                
                MDeltaAngularVelocity.X += angularComponent.X * (impulseMagnitude * MAngularFactor.X);
                MDeltaAngularVelocity.Y += angularComponent.Y  * (impulseMagnitude * MAngularFactor.Y);
                MDeltaAngularVelocity.Z += angularComponent.Z *(impulseMagnitude * MAngularFactor.Z);

                //m_deltaAngularVelocity += angularComponent*(impulseMagnitude*m_angularFactor);


                MathUtil.SanityCheckVector(ref MDeltaLinearVelocity);
                MathUtil.SanityCheckVector(ref MDeltaAngularVelocity);
            }

	    }


        public void InternalApplyPushImpulse(IndexedVector3 linearComponent, IndexedVector3 angularComponent, float impulseMagnitude)
        {
            InternalApplyPushImpulse(ref linearComponent, ref angularComponent, impulseMagnitude);
        }

        public void InternalApplyPushImpulse(ref IndexedVector3 linearComponent, ref IndexedVector3 angularComponent,float impulseMagnitude)
	    {
		    if (_mInverseMass != 0f)
		    {
			    MPushVelocity += linearComponent*impulseMagnitude;
			    MTurnVelocity += angularComponent*(impulseMagnitude*MAngularFactor);
		    }
	    }
    	
	    public void	InternalWritebackVelocity()
	    {
		    if (_mInverseMass != 0f)
		    {
			    SetLinearVelocity(GetLinearVelocity()+ MDeltaLinearVelocity);
			    SetAngularVelocity(GetAngularVelocity()+MDeltaAngularVelocity);
			    //m_deltaLinearVelocity = IndexedVector3.Zero;
                //m_deltaAngularVelocity = IndexedVector3.Zero;
			    //m_originalBody->setCompanionId(-1);
		    }
	    }


        

        public void InternalWritebackVelocity(float timeStep)
        {
#if DEBUG        
			if (BulletGlobals.g_streamWriter != null && BulletGlobals.debugRigidBody)
			{
				BulletGlobals.g_streamWriter.WriteLine(String.Format("[{0}] internalWritebackVelocity ",(String)m_userObjectPointer));
				MathUtil.PrintVector3(BulletGlobals.g_streamWriter,"LinearVelocity",GetLinearVelocity());
				MathUtil.PrintVector3(BulletGlobals.g_streamWriter, "DeltaLinearVelocity", MDeltaLinearVelocity);
				MathUtil.PrintVector3(BulletGlobals.g_streamWriter, "AngularVelocity", GetAngularVelocity());
				MathUtil.PrintVector3(BulletGlobals.g_streamWriter, "DeltaAngularVelocity", MDeltaAngularVelocity);
			}
#endif
	        if (_mInverseMass != 0f)
	        {
		        SetLinearVelocity(GetLinearVelocity()+ MDeltaLinearVelocity);
		        SetAngularVelocity(GetAngularVelocity()+MDeltaAngularVelocity);
        		
		        //correct the position/orientation based on push/turn recovery
		        IndexedMatrix newTransform;
		        TransformUtil.IntegrateTransform(GetWorldTransform(),MPushVelocity,MTurnVelocity,timeStep,out newTransform);
		        SetWorldTransform(ref newTransform);
		        //m_originalBody->setCompanionId(-1);
	        }
#if DEBUG
			if (BulletGlobals.g_streamWriter != null && BulletGlobals.debugRigidBody)
            {
                BulletGlobals.g_streamWriter.WriteLine("post integrate transform.");
                MathUtil.PrintMatrix(BulletGlobals.g_streamWriter, GetWorldTransform());
            }
#endif
			//m_deltaLinearVelocity = IndexedVector3.Zero;
			//m_deltaAngularVelocity = IndexedVector3.Zero;
			//m_pushVelocity = IndexedVector3.Zero;
			//m_turnVelocity = IndexedVector3.Zero;
        }
    }


	///The btRigidBodyConstructionInfo structure provides information to create a rigid body. Setting mass to zero creates a fixed (non-dynamic) rigid body.
	///For dynamic objects, you can use the collision shape to approximate the local inertia tensor, otherwise use the zero vector (default argument)
	///You can use the motion state to synchronize the world transform between physics and graphics objects. 
	///And if the motion state is provided, the rigid body will initialize its initial world transform from the motion state,
	///m_startWorldTransform is only used when you don't provide a motion state.
	public class RigidBodyConstructionInfo
	{
		public float m_mass;

		///When a motionState is provided, the rigid body will initialize its world transform from the motion state
		///In this case, m_startWorldTransform is ignored.
		public IMotionState		m_motionState;
		public IndexedMatrix	m_startWorldTransform;

		public CollisionShape	m_collisionShape;
		public IndexedVector3			m_localInertia;
		public float			m_linearDamping;
		public float			m_angularDamping;

		///best simulation results when friction is non-zero
		public float			m_friction;
		///best simulation results using zero restitution.
		public float			m_restitution;

		public float			m_linearSleepingThreshold;
		public float			m_angularSleepingThreshold;

		//Additional damping can help avoiding lowpass jitter motion, help stability for ragdolls etc.
		//Such damping is undesirable, so once the overall simulation quality of the rigid body dynamics system has improved, this should become obsolete
		public bool				m_additionalDamping;
		public float			m_additionalDampingFactor;
		public float			m_additionalLinearDampingThresholdSqr;
		public float			m_additionalAngularDampingThresholdSqr;
		public float			m_additionalAngularDampingFactor;

        public RigidBodyConstructionInfo(float mass, IMotionState motionState, CollisionShape collisionShape): this(mass,motionState,collisionShape,new IndexedVector3(0))
        {

        }
		public RigidBodyConstructionInfo(float mass, IMotionState motionState, CollisionShape collisionShape, IndexedVector3 localInertia)
        {
    		m_mass = mass;
			m_motionState =motionState;
			m_collisionShape = collisionShape;
			m_localInertia = localInertia;
			m_linearDamping = 0f;
			m_angularDamping = 0f;
			m_friction = 0.5f;
			m_restitution = 0f;
			m_linearSleepingThreshold = 0.8f;
			m_angularSleepingThreshold = 1f;
			m_additionalDamping = false;
			m_additionalDampingFactor = 0.005f;
			m_additionalLinearDampingThresholdSqr = 0.01f;
			m_additionalAngularDampingThresholdSqr = 0.01f;
			m_additionalAngularDampingFactor = 0.01f;
            m_startWorldTransform = IndexedMatrix.Identity;
		}
	}

    [Flags]
    public enum RigidBodyFlags
    {
        BT_NONE = 0,
        BT_DISABLE_WORLD_GRAVITY = 1
    }



}
