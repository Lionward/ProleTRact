import React, { useState, useEffect } from 'react';
import './TandemRepeatAnimation.css';

const TandemRepeatAnimation: React.FC = () => {
  const [animationStep, setAnimationStep] = useState(0);
  const [isPlaying, setIsPlaying] = useState(true);

  const steps = [
    { title: 'Normal Sequence', repeats: 3, description: 'A normal sequence with 3 CAG repeats' },
    { title: 'Expanded Sequence', repeats: 8, description: 'An expanded sequence with 8 CAG repeats' },
    { title: 'Pathogenic Expansion', repeats: 15, description: 'A pathogenic expansion with 15+ CAG repeats' }
  ];

  useEffect(() => {
    if (!isPlaying) return;

    const interval = setInterval(() => {
      setAnimationStep((prev) => (prev + 1) % steps.length);
    }, 3000);

    return () => clearInterval(interval);
  }, [isPlaying, steps.length]);

  const currentStep = steps[animationStep];
  const motif = 'CAG';
  const motifColor = '#667eea';
  const interruptionColor = '#4b5563';

  // Create sequence with repeats
  const createSequence = (repeatCount: number) => {
    const sequence = [];
    for (let i = 0; i < repeatCount; i++) {
      sequence.push({ type: 'motif', text: motif, index: i });
    }
    return sequence;
  };

  const sequence = createSequence(currentStep.repeats);

  return (
    <div className="tandem-repeat-animation">
      <div className="animation-controls">
        <button 
          className="play-pause-btn"
          onClick={() => setIsPlaying(!isPlaying)}
          aria-label={isPlaying ? 'Pause animation' : 'Play animation'}
        >
          {isPlaying ? '⏸️' : '▶️'}
        </button>
        <div className="step-indicators">
          {steps.map((_, index) => (
            <button
              key={index}
              className={`step-dot ${index === animationStep ? 'active' : ''}`}
              onClick={() => {
                setAnimationStep(index);
                setIsPlaying(false);
              }}
              aria-label={`Go to step ${index + 1}`}
            />
          ))}
        </div>
      </div>

      <div className="animation-container">
        <div className="dna-strand">
          <div className="dna-backbone left"></div>
          <div className="dna-sequence">
            <div className="sequence-label">Tandem Repeat Region</div>
            <div className="sequence-display">
              {sequence.map((segment, idx) => (
                <div
                  key={idx}
                  className={`sequence-segment motif-segment`}
                  style={{
                    backgroundColor: motifColor,
                    animationDelay: `${idx * 0.1}s`
                  }}
                >
                  <span className="base-text">{segment.text}</span>
                </div>
              ))}
            </div>
            <div className="repeat-count">
              <span className="count-label">Repeat Count:</span>
              <span className="count-value">{currentStep.repeats}</span>
            </div>
          </div>
          <div className="dna-backbone right"></div>
        </div>

        <div className="animation-info">
          <h3 className="step-title">{currentStep.title}</h3>
          <p className="step-description">{currentStep.description}</p>
        </div>

        <div className="expansion-visualization">
          <div className="expansion-bar">
            <div 
              className="expansion-fill"
              style={{ 
                width: `${(currentStep.repeats / 15) * 100}%`,
                backgroundColor: currentStep.repeats < 5 
                  ? '#10b981' 
                  : currentStep.repeats < 10 
                  ? '#f59e0b' 
                  : '#ef4444'
              }}
            >
              <span className="expansion-label">
                {currentStep.repeats < 5 
                  ? 'Normal' 
                  : currentStep.repeats < 10 
                  ? 'Intermediate' 
                  : 'Pathogenic'}
              </span>
            </div>
          </div>
          <div className="expansion-markers">
            <span className="marker">Normal (&lt;5)</span>
            <span className="marker">Intermediate (5-10)</span>
            <span className="marker">Pathogenic (&gt;10)</span>
          </div>
        </div>
      </div>

      <div className="explanation-box">
        <div className="explanation-icon">💡</div>
        <div className="explanation-content">
          <h4>What are Tandem Repeats?</h4>
          <p>
            Tandem repeats are sequences of DNA that are repeated multiple times in a row. 
            When these repeats expand beyond normal ranges, they can cause genetic disorders. 
            ProleTRact helps you visualize and analyze these expansions to understand their impact.
          </p>
        </div>
      </div>
    </div>
  );
};

export default TandemRepeatAnimation;

