import React, { useState, useEffect } from 'react';
import './AnnotationPanel.css';

export interface Annotation {
  id: string;
  region: string;
  comment: string;
  tags: string[];
  color: string;
  timestamp: number;
}

interface AnnotationPanelProps {
  region: string;
  onAnnotationChange?: (annotations: Annotation[]) => void;
}

const STORAGE_KEY = 'proletract_annotations';
const TAG_COLORS = ['#667eea', '#764ba2', '#f093fb', '#f5576c', '#4fd1c7', '#68d391', '#f6e05e', '#f6ad55'];

const AnnotationPanel: React.FC<AnnotationPanelProps> = ({ region, onAnnotationChange }) => {
  const [annotations, setAnnotations] = useState<Annotation[]>([]);
  const [isOpen, setIsOpen] = useState(false);
  const [comment, setComment] = useState('');
  const [selectedTags, setSelectedTags] = useState<string[]>([]);
  const [tagInput, setTagInput] = useState('');
  const [availableTags, setAvailableTags] = useState<string[]>([]);

  // Load annotations from localStorage
  useEffect(() => {
    try {
      const stored = localStorage.getItem(STORAGE_KEY);
      if (stored) {
        const allAnnotations: Annotation[] = JSON.parse(stored);
        const regionAnnotations = allAnnotations.filter(a => a.region === region);
        setAnnotations(regionAnnotations);
        
        // Extract all unique tags
        const allTags = new Set<string>();
        allAnnotations.forEach(a => a.tags.forEach(tag => allTags.add(tag)));
        setAvailableTags(Array.from(allTags));
      }
    } catch (error) {
      console.error('Error loading annotations:', error);
    }
  }, [region]);

  // Save annotations to localStorage
  useEffect(() => {
    try {
      const stored = localStorage.getItem(STORAGE_KEY);
      const allAnnotations: Annotation[] = stored ? JSON.parse(stored) : [];
      const otherAnnotations = allAnnotations.filter(a => a.region !== region);
      const updated = [...otherAnnotations, ...annotations];
      localStorage.setItem(STORAGE_KEY, JSON.stringify(updated));
      
      if (onAnnotationChange) {
        onAnnotationChange(annotations);
      }
    } catch (error) {
      console.error('Error saving annotations:', error);
    }
  }, [annotations, region, onAnnotationChange]);

  const handleAddAnnotation = () => {
    if (!comment.trim() && selectedTags.length === 0) return;

    const newAnnotation: Annotation = {
      id: `annotation_${Date.now()}_${Math.random().toString(36).substr(2, 9)}`,
      region,
      comment: comment.trim(),
      tags: [...selectedTags],
      color: TAG_COLORS[annotations.length % TAG_COLORS.length],
      timestamp: Date.now(),
    };

    setAnnotations(prev => [...prev, newAnnotation]);
    setComment('');
    setSelectedTags([]);
    setIsOpen(false);
  };

  const handleDeleteAnnotation = (id: string) => {
    setAnnotations(prev => prev.filter(a => a.id !== id));
  };

  const handleAddTag = () => {
    const tag = tagInput.trim().toLowerCase();
    if (tag && !selectedTags.includes(tag) && !availableTags.includes(tag)) {
      setAvailableTags(prev => [...prev, tag]);
    }
    if (tag && !selectedTags.includes(tag)) {
      setSelectedTags(prev => [...prev, tag]);
    }
    setTagInput('');
  };

  const toggleTag = (tag: string) => {
    setSelectedTags(prev =>
      prev.includes(tag) ? prev.filter(t => t !== tag) : [...prev, tag]
    );
  };

  const regionAnnotations = annotations.filter(a => a.region === region);

  return (
    <div className="annotation-panel">
      <button
        className="annotation-toggle-button"
        onClick={() => setIsOpen(!isOpen)}
        title="Add annotation"
      >
        <span className="annotation-icon">📝</span>
        <span className="annotation-label">Annotations</span>
        {regionAnnotations.length > 0 && (
          <span className="annotation-count">{regionAnnotations.length}</span>
        )}
      </button>

      {isOpen && (
        <div className="annotation-dialog">
          <div className="annotation-dialog-header">
            <h3>Add Annotation for {region}</h3>
            <button
              className="annotation-close-button"
              onClick={() => setIsOpen(false)}
            >
              ✕
            </button>
          </div>

          <div className="annotation-dialog-content">
            <div className="annotation-comment-section">
              <label>Comment:</label>
              <textarea
                className="annotation-comment-input"
                value={comment}
                onChange={(e) => setComment(e.target.value)}
                placeholder="Add a comment about this region..."
                rows={4}
              />
            </div>

            <div className="annotation-tags-section">
              <label>Tags:</label>
              <div className="annotation-tag-input-group">
                <input
                  type="text"
                  className="annotation-tag-input"
                  value={tagInput}
                  onChange={(e) => setTagInput(e.target.value)}
                  onKeyPress={(e) => {
                    if (e.key === 'Enter') {
                      e.preventDefault();
                      handleAddTag();
                    }
                  }}
                  placeholder="Add tag and press Enter"
                />
                <button
                  className="annotation-tag-add-button"
                  onClick={handleAddTag}
                >
                  Add
                </button>
              </div>

              {availableTags.length > 0 && (
                <div className="annotation-tags-list">
                  {availableTags.map(tag => (
                    <button
                      key={tag}
                      className={`annotation-tag ${selectedTags.includes(tag) ? 'selected' : ''}`}
                      onClick={() => toggleTag(tag)}
                    >
                      {tag}
                    </button>
                  ))}
                </div>
              )}
            </div>

            <button
              className="annotation-save-button"
              onClick={handleAddAnnotation}
              disabled={!comment.trim() && selectedTags.length === 0}
            >
              Save Annotation
            </button>
          </div>
        </div>
      )}

      {regionAnnotations.length > 0 && (
        <div className="annotation-list">
          {regionAnnotations.map(annotation => (
            <div key={annotation.id} className="annotation-item">
              <div className="annotation-item-header">
                <span className="annotation-item-color" style={{ backgroundColor: annotation.color }} />
                <span className="annotation-item-time">
                  {new Date(annotation.timestamp).toLocaleString()}
                </span>
                <button
                  className="annotation-item-delete"
                  onClick={() => handleDeleteAnnotation(annotation.id)}
                  title="Delete annotation"
                >
                  🗑️
                </button>
              </div>
              {annotation.comment && (
                <div className="annotation-item-comment">{annotation.comment}</div>
              )}
              {annotation.tags.length > 0 && (
                <div className="annotation-item-tags">
                  {annotation.tags.map(tag => (
                    <span key={tag} className="annotation-item-tag">
                      {tag}
                    </span>
                  ))}
                </div>
              )}
            </div>
          ))}
        </div>
      )}
    </div>
  );
};

export default AnnotationPanel;




